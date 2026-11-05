//! Local source dispatch and no-clobber publication for selection indexes.

use std::{
    fs::{self, File},
    io::{self, BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
    sync::Arc,
};

use tempfile::NamedTempFile;

use super::{ErrorRecordingReader, SaveIndex, SelectionIndex, construction::construct_reference};
use crate::error::{IoError, Result, UnavailableInputError, UsageError, constructors};

pub(super) struct Source<'a> {
    path: &'a Path,
    role: &'static str,
    reader: BufReader<ErrorRecordingReader<File>>,
    reference: bool,
    save: Option<(PathBuf, NamedTempFile)>,
}

impl<'a> Source<'a> {
    pub(super) fn open(
        path: &'a Path,
        role: &'static str,
        save: Option<SaveIndex<'_>>,
        k: u8,
        w: u8,
    ) -> Result<Self> {
        let file = File::open(path).map_err(|source| UnavailableInputError::SelectionIndex {
            role,
            path: path.to_path_buf(),
            source,
        })?;
        let mut reader = BufReader::with_capacity(
            1 << 20,
            ErrorRecordingReader {
                inner: file,
                error_kind: None,
            },
        );
        let prefix = reader
            .fill_buf()
            .map_err(|source| IoError::SelectionIndex {
                role,
                path: path.to_path_buf(),
                source: Box::new(source),
            })?;
        // Match needletail's compression probes. A recognized prefix commits to
        // reference parsing; a decode failure is never retried as an index.
        let reference = matches!(prefix.first(), Some(b'>' | b'@'))
            || [b"\x1f\x8b", b"BZ", b"\xfd7", b"\x28\xb5"]
                .iter()
                .any(|magic| prefix.starts_with(*magic));
        if save.is_some() && !reference {
            return Err(UsageError::SaveRequiresReference { role }.into());
        }
        let save = save
            .map(|requested| {
                let destination = match requested {
                    SaveIndex::BesideReference => {
                        let mut name = path.as_os_str().to_os_string();
                        name.push(format!(".k{k}.w{w}.idx"));
                        PathBuf::from(name)
                    }
                    SaveIndex::At(path) => path.to_path_buf(),
                };
                ensure_absent(&destination)?;
                let parent = parent_dir(&destination);
                // Fix the destination's directory identity before construction.
                let parent =
                    fs::canonicalize(parent).map_err(|source| save_error(&destination, source))?;
                let destination = parent.join(destination.file_name().ok_or_else(|| {
                    save_error(
                        &destination,
                        io::Error::new(
                            io::ErrorKind::InvalidInput,
                            "index destination needs a filename",
                        ),
                    )
                })?);
                let staged = tempfile::Builder::new()
                    .prefix(".nuclease-index-")
                    .tempfile_in(parent)
                    .map_err(|source| save_error(&destination, source))?;
                Ok::<_, crate::error::RunError>((destination, staged))
            })
            .transpose()?;
        Ok(Self {
            path,
            role,
            reader,
            reference,
            save,
        })
    }

    pub(super) fn is_reference(&self) -> bool {
        self.reference
    }

    pub(super) fn save_path(&self) -> Option<&Path> {
        self.save.as_ref().map(|(path, _)| path.as_path())
    }

    pub(super) fn prepare(mut self, k: u8, w: u8, min_hits: usize) -> Result<SelectionIndex> {
        let index = if self.reference {
            construct_reference(self.reader, self.role, self.path, k, w)?
        } else {
            deacon::load_index(&mut self.reader).map_err(|source| {
                constructors::index_load_error(
                    self.role,
                    self.path,
                    self.reader.get_ref().error_kind,
                    source.into_boxed_dyn_error(),
                )
            })?
        };
        if let Some((path, mut staged)) = self.save {
            {
                let mut writer = BufWriter::new(&mut staged);
                deacon::write_index(&index, &mut writer).map_err(|source| IoError::SaveIndex {
                    path: path.clone(),
                    source: source.into_boxed_dyn_error(),
                })?;
                writer.flush().map_err(|source| save_error(&path, source))?;
            }
            // Only complete bytes become visible at the destination. Persisting
            // never replaces a file created by another process during the build.
            staged
                .persist_noclobber(&path)
                .map_err(|error| save_error(&path, error.error))?;
        }
        Ok(SelectionIndex {
            index: Arc::new(index),
            min_hits,
        })
    }
}

pub(super) fn reject_same_sources(enrich: Option<&Path>, deplete: Option<&Path>) -> Result<()> {
    if let (Some(enrich), Some(deplete)) = (enrich, deplete) {
        let identity = |path: &Path, role| {
            same_file::Handle::from_path(path).map_err(|source| {
                UnavailableInputError::SelectionIndex {
                    role,
                    path: path.to_path_buf(),
                    source,
                }
            })
        };
        if identity(enrich, "enrichment")? == identity(deplete, "depletion")? {
            return Err(UsageError::SameSelectionSource {
                enrich: enrich.to_path_buf(),
                deplete: deplete.to_path_buf(),
            }
            .into());
        }
    }
    Ok(())
}

fn parent_dir(path: &Path) -> &Path {
    path.parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or(Path::new("."))
}

fn save_error(path: &Path, source: io::Error) -> crate::error::RunError {
    IoError::SaveIndex {
        path: path.to_path_buf(),
        source: Box::new(source),
    }
    .into()
}

fn ensure_absent(path: &Path) -> Result<()> {
    match fs::symlink_metadata(path) {
        Ok(_) => Err(UsageError::IndexDestinationExists {
            path: path.to_path_buf(),
        }
        .into()),
        Err(error) if error.kind() == io::ErrorKind::NotFound => Ok(()),
        Err(error) => Err(save_error(path, error)),
    }
}

// Resolve an output even when its leaf is a dangling symlink to a future index.
// A missing parent cannot alias a save, whose parent must already exist.
fn output_target(path: &Path) -> io::Result<Option<PathBuf>> {
    let mut path = path.to_path_buf();
    for _ in 0..40 {
        match fs::symlink_metadata(&path) {
            Ok(metadata) if metadata.file_type().is_symlink() => {
                let target = fs::read_link(&path)?;
                path = parent_dir(&path).join(target);
            }
            Ok(_) => return fs::canonicalize(path).map(Some),
            Err(error) if error.kind() == io::ErrorKind::NotFound => {
                return match fs::canonicalize(parent_dir(&path)) {
                    Ok(parent) => Ok(path.file_name().map(|name| parent.join(name))),
                    Err(error) if error.kind() == io::ErrorKind::NotFound => Ok(None),
                    Err(error) => Err(error),
                };
            }
            Err(error) => return Err(error),
        }
    }
    Err(io::Error::new(
        io::ErrorKind::InvalidInput,
        "too many output symlink levels",
    ))
}

pub(super) fn validate_destinations(
    saves: &[Option<&Path>],
    outputs: &[Option<&Path>],
) -> Result<()> {
    for (position, path) in saves
        .iter()
        .enumerate()
        .filter_map(|(i, path)| path.map(|path| (i, path)))
    {
        for other in saves[..position].iter().chain(outputs).flatten() {
            let target = output_target(other).map_err(|source| save_error(path, source))?;
            if target
                .as_deref()
                .map(|target| paths_collide(path, target))
                .transpose()
                .map_err(|source| save_error(path, source))?
                .unwrap_or(false)
            {
                return Err(UsageError::IndexOutputCollision {
                    path: path.to_path_buf(),
                    other: other.to_path_buf(),
                }
                .into());
            }
        }
    }
    Ok(())
}

fn paths_collide(save: &Path, other: &Path) -> io::Result<bool> {
    if save == other {
        return Ok(true);
    }
    if !same_file::is_same_file(parent_dir(save), parent_dir(other))? {
        return Ok(false);
    }
    let Some(other_name) = other.file_name() else {
        return Ok(false);
    };
    // Absent leaf names cannot be compared by file identity yet. Let this
    // filesystem apply its case/Unicode rules in an isolated temporary directory,
    // rather than guessing those rules or exposing an incomplete destination.
    let probe = tempfile::Builder::new()
        .prefix(".nuclease-path-check-")
        .tempdir_in(parent_dir(save))?;
    let original = probe
        .path()
        .join(save.file_name().expect("validated save filename"));
    let alternate = probe.path().join(other_name);
    File::create_new(&original)?;
    match same_file::is_same_file(original, alternate) {
        Ok(collides) => Ok(collides),
        Err(error) if error.kind() == io::ErrorKind::NotFound => Ok(false),
        Err(error) => Err(error),
    }
}

pub(super) fn validate_published_outputs(
    saves: &[Option<&Path>],
    outputs: &[Option<&Path>],
) -> Result<()> {
    for path in saves.iter().flatten() {
        for other in outputs.iter().flatten() {
            match same_file::is_same_file(path, other) {
                Ok(true) => {
                    return Err(UsageError::IndexOutputCollision {
                        path: path.to_path_buf(),
                        other: other.to_path_buf(),
                    }
                    .into());
                }
                Ok(false) => {}
                Err(error) if error.kind() == io::ErrorKind::NotFound => {}
                Err(error) => return Err(save_error(path, error)),
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::fs;

    use super::Source;
    use crate::filter::minimizer::SaveIndex;

    #[test]
    fn publication_does_not_replace_a_destination_created_after_preflight() {
        let dir = tempfile::tempdir().expect("tempdir");
        let reference = dir.path().join("reference.fa");
        let destination = dir.path().join("saved.idx");
        fs::write(&reference, b">r\nACGTGCTAACGT\n").expect("reference");
        let source = Source::open(
            &reference,
            "enrichment",
            Some(SaveIndex::At(&destination)),
            3,
            1,
        )
        .expect("preflight");
        fs::write(&destination, b"other process").expect("competing destination");
        assert!(source.prepare(3, 1, 1).is_err());
        assert_eq!(
            fs::read(&destination).expect("destination"),
            b"other process"
        );
        assert_eq!(fs::read_dir(dir.path()).expect("entries").count(), 2);
    }
}
