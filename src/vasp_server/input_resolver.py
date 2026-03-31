import shutil
from pathlib import Path
from typing import Any


class UpstreamInputResolver:
    def materialize_inputs(self, task_type: str, upstream_artifacts: list[dict], work_dir: Path) -> dict[str, str]:
        if task_type == "scf_calculation":
            return self.resolve_for_scf(upstream_artifacts, work_dir)
        if task_type == "dos_calculation":
            return self.resolve_for_dos(upstream_artifacts, work_dir)
        raise ValueError(f"unsupported task type for upstream resolution: {task_type}")

    def resolve_for_scf(self, artifacts: list[dict], work_dir: Path) -> dict[str, str]:
        source = self._find_first(artifacts, ["contcar", "poscar"])
        if source is None:
            raise FileNotFoundError("上游产物中缺少 CONTCAR 或 POSCAR")

        destination = work_dir / "POSCAR"
        self._materialize_artifact(source, destination)
        return {"POSCAR": str(destination)}

    def resolve_for_dos(self, artifacts: list[dict], work_dir: Path) -> dict[str, str]:
        copied_files: dict[str, str] = {}
        mappings = [
            ("POSCAR", ["poscar", "contcar"]),
            ("POTCAR", ["potcar"]),
            ("CHG", ["chg"]),
            ("CHGCAR", ["chgcar"]),
            ("WAVECAR", ["wavecar"]),
        ]

        for target_name, artifact_types in mappings:
            source = self._find_first(artifacts, artifact_types)
            if source is None:
                if target_name in {"POSCAR", "POTCAR", "CHGCAR"}:
                    raise FileNotFoundError(f"上游产物中缺少关键文件 {target_name}")
                continue
            destination = work_dir / target_name
            self._materialize_artifact(source, destination)
            copied_files[target_name] = str(destination)

        return copied_files

    def _find_first(self, artifacts: list[dict], artifact_types: list[str]) -> dict[str, Any] | None:
        normalized = {artifact_type.lower() for artifact_type in artifact_types}
        for artifact in artifacts:
            artifact_type = str(artifact.get("artifact_type", "")).lower()
            filename = str(artifact.get("filename", "")).lower()
            if artifact_type in normalized or filename in normalized:
                return artifact
        return None

    def _materialize_artifact(self, artifact: dict[str, Any], destination: Path) -> None:
        source_path = artifact.get("local_path") or artifact.get("source_path") or artifact.get("storage_key")
        if not source_path:
            raise FileNotFoundError(f"产物缺少可用源路径: {artifact}")

        src = Path(str(source_path))
        if not src.exists():
            raise FileNotFoundError(f"产物源文件不存在: {src}")

        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(str(src), str(destination))
