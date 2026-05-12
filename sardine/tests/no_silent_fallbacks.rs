//! Forbidden-pattern guard.
//!
//! Runs `scripts/check_no_silent_fallbacks.sh` against `src/` and asserts
//! that no unannotated silent-fallback patterns exist.  See
//! [`AGENTS.md`](../../../AGENTS.md) for the policy and the SAFETY-OK
//! annotation convention.

use std::path::PathBuf;
use std::process::Command;

/// Locate the repository root by walking upward until we find `AGENTS.md`.
fn repo_root() -> PathBuf {
    let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    loop {
        if p.join("AGENTS.md").is_file() {
            return p;
        }
        if !p.pop() {
            panic!("could not locate repository root from CARGO_MANIFEST_DIR");
        }
    }
}

#[test]
fn no_unannotated_silent_fallbacks_in_src() {
    let root = repo_root();
    let script = root.join("scripts/check_no_silent_fallbacks.sh");
    assert!(
        script.is_file(),
        "guard script missing at {}",
        script.display()
    );

    let output = Command::new("bash")
        .arg(&script)
        .arg("sardine/src")
        .current_dir(&root)
        .output()
        .expect("failed to invoke check_no_silent_fallbacks.sh");

    if !output.status.success() {
        let stdout = String::from_utf8_lossy(&output.stdout);
        let stderr = String::from_utf8_lossy(&output.stderr);
        panic!(
            "forbidden-pattern guard failed.\n\n--- stdout ---\n{}\n--- stderr ---\n{}",
            stdout, stderr
        );
    }
}
