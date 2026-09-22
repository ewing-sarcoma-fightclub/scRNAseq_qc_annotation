"""Offline orchestration checks; these do not validate biological results."""

import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


def executable(path, body):
    path.write_text(f"#!{sys.executable}\n" + body)
    path.chmod(0o755)
    return path


class ConfigurationTests(unittest.TestCase):
    def test_config_reaches_analysis_processes_and_disables_resume(self):
        with tempfile.TemporaryDirectory(prefix="scrna config ") as tmp:
            work = Path(tmp)
            inputs = work / "inputs"
            (inputs / "sample1").mkdir(parents=True)
            log = work / "calls.jsonl"
            fake_r = executable(work / "Rscript", """
import json, os, sys
with open(os.environ['CHECK_LOG'], 'a') as handle:
    handle.write(json.dumps({'script': os.path.basename(sys.argv[1]),
                            'env': dict(os.environ)}) + '\\n')
""")
            config = work / "config.env"
            config.write_text((ROOT / "env/config.env").read_text() + "\n" +
                              f"R_BIN={shlex.quote(str(fake_r))}\n" +
                              "RESUME=false\nANNOTATION_DE_LOGFC=0.73\n")
            env = dict(os.environ, CHECK_LOG=str(log))
            keys = ["AUCELL_SLOT", "DOUBLET_FILTER_MODE", "ANNOTATION_TRIAGE_ENABLE",
                    "ANNOTATION_DE_LOGFC", "RESUME"]
            for key in keys:
                env.pop(key, None)
            for entrypoint in ("run_all.sh", "pipeline_QC_after_cellranger.sh"):
                with self.subTest(entrypoint=entrypoint):
                    out = work / entrypoint / "qc"
                    previous = out / "emptydrops/sample1/sample1_emptydrops_results.csv"
                    previous.parent.mkdir(parents=True)
                    previous.write_text("old result\n")
                    log.write_text("")
                    if entrypoint == "run_all.sh":
                        args = ["--skip-cellranger", "--cellranger-root", str(inputs),
                                "--qc-out", str(out)]
                    else:
                        args = ["--root", str(inputs), "--out", str(out)]
                    result = subprocess.run(
                        ["bash", str(ROOT / "bin" / entrypoint), *args,
                         "--config", str(config)], env=env, text=True, capture_output=True)
                    self.assertEqual(result.returncode, 0, result.stderr)
                    calls = [json.loads(line) for line in log.read_text().splitlines()]
                    self.assertIn("EmptyDrops_per_sample.R", [x["script"] for x in calls])
                    self.assertIn("Seurat_merge_annotate_integrate.R", [x["script"] for x in calls])
                    for call in calls:
                        for key, expected in {"AUCELL_SLOT": "counts",
                                              "DOUBLET_FILTER_MODE": "union",
                                              "ANNOTATION_TRIAGE_ENABLE": "true",
                                              "ANNOTATION_DE_LOGFC": "0.73",
                                              "RESUME": "false"}.items():
                            self.assertEqual(call["env"].get(key), expected, (call["script"], key))

    def test_setup_in_fresh_environment_without_conda_envs_path(self):
        with tempfile.TemporaryDirectory(prefix="scrna setup ") as tmp:
            work = Path(tmp)
            fake_mamba = executable(work / "micromamba", """
import os, pathlib, sys
root = pathlib.Path(os.environ['CHECK_ENVS'])
if sys.argv[1:3] == ['env', 'list']:
    for path in sorted(root.glob('*')):
        print(path.name, str(path))
elif sys.argv[1:3] == ['env', 'create']:
    name = sys.argv[sys.argv.index('-n') + 1]
    (root / name / 'bin').mkdir(parents=True)
else:
    raise SystemExit('Unexpected fake-mamba command')
""")
            env = dict(os.environ, CHECK_ENVS=str(work / "envs"))
            env.pop("CONDA_ENVS_PATH", None)
            env.pop("MAMBA_ROOT_PREFIX", None)
            config = work / "generated.env"
            result = subprocess.run(
                ["bash", str(ROOT / "bin/setup_envs.sh"), "--mamba-bin", str(fake_mamba),
                 "--skip-r", "--config-out", str(config)],
                env=env, text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertNotIn("unbound variable", result.stderr)
            self.assertIn(str(work / "envs/ewing-scrna-py/bin/python"), config.read_text())


if __name__ == "__main__":
    unittest.main()
