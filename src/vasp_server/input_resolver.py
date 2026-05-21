import shutil
from pathlib import Path
from typing import Any
import requests


class UpstreamInputResolver:
    def materialize_inputs(self, task_type: str, upstream_artifacts: list[dict], work_dir: Path) -> dict[str, str]:
        if task_type == "scf_calculation":
            return self.resolve_for_scf(upstream_artifacts, work_dir)
        if task_type == "dos_calculation":
            return self.resolve_for_dos(upstream_artifacts, work_dir)
        if task_type == "band_structure":
            return self.resolve_for_band_structure(upstream_artifacts, work_dir)
        if task_type == "md_calculation":
            return self.resolve_for_md(upstream_artifacts, work_dir)
        if task_type == "phonon_calculation":
            return self.resolve_for_phonon(upstream_artifacts, work_dir)
        if task_type == "custom_calculation":
            return self.resolve_for_custom(upstream_artifacts, work_dir)
        raise ValueError(f"unsupported task type for upstream resolution: {task_type}")

    def resolve_for_scf(self, artifacts: list[dict], work_dir: Path) -> dict[str, str]:
        source = self._find_first(artifacts, ["contcar", "poscar"])
        if source is None:
            raise FileNotFoundError("上游产物中缺少 CONTCAR 或 POSCAR")

        destination = work_dir / "POSCAR"
        self._materialize_artifact(source, destination, require_nonempty=True)
        return {"POSCAR": str(destination)}

    def resolve_single_structure(self, artifacts: list[dict], work_dir: Path, dest_name: str = "POSCAR") -> dict[str, str]:
        source = self._find_first(artifacts, ["contcar", "poscar"])
        if source is None:
            raise FileNotFoundError("上游产物中缺少 CONTCAR 或 POSCAR")

        destination = work_dir / dest_name
        self._materialize_artifact(source, destination, require_nonempty=True)
        return {dest_name: str(destination)}

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
            self._materialize_artifact(
                source,
                destination,
                require_nonempty=target_name in {"POSCAR", "POTCAR", "CHGCAR"},
                empty_error=(
                    "上游产物中的 CHGCAR 为空，未生成有效的 CHGCAR"
                    if target_name == "CHGCAR"
                    else f"上游产物中的 {target_name} 为空"
                ),
            )
            copied_files[target_name] = str(destination)

        return copied_files

    def resolve_for_band_structure(self, artifacts: list[dict], work_dir: Path) -> dict[str, str]:
        copied_files: dict[str, str] = {}
        mappings = [
            ("POSCAR", ["contcar", "poscar"]),
            ("POTCAR", ["potcar"]),
            ("CHG", ["chg"]),
            ("CHGCAR", ["chgcar"]),
        ]
        required = {"POSCAR", "POTCAR", "CHGCAR"}
        for target_name, artifact_types in mappings:
            source = self._find_first(artifacts, artifact_types)
            if source is None:
                if target_name in required:
                    raise FileNotFoundError(f"上游产物中缺少关键文件 {target_name}")
                continue
            destination = work_dir / target_name
            self._materialize_artifact(
                source,
                destination,
                require_nonempty=target_name in required,
                empty_error=(
                    "上游产物中的 CHGCAR 为空，未生成有效的 CHGCAR"
                    if target_name == "CHGCAR"
                    else f"上游产物中的 {target_name} 为空"
                ),
            )
            copied_files[target_name] = str(destination)
        return copied_files

    def resolve_for_md(self, artifacts: list[dict], work_dir: Path) -> dict[str, str]:
        copied_files: dict[str, str] = {}
        mappings = [
            ("POSCAR", ["contcar", "poscar"]),
            ("POTCAR", ["potcar"]),
        ]
        for target_name, artifact_types in mappings:
            source = self._find_first(artifacts, artifact_types)
            if source is None:
                raise FileNotFoundError(f"上游产物中缺少关键文件 {target_name}")
            destination = work_dir / target_name
            self._materialize_artifact(source, destination, require_nonempty=True)
            copied_files[target_name] = str(destination)
        return copied_files

    def resolve_for_phonon(self, artifacts: list[dict], work_dir: Path) -> dict[str, str]:
        return self.resolve_for_md(artifacts, work_dir)

    def resolve_for_custom(self, artifacts: list[dict], work_dir: Path) -> dict[str, str]:
        return self.resolve_single_structure(artifacts, work_dir, dest_name="POSCAR")

    def _find_first(self, artifacts: list[dict], artifact_types: list[str]) -> dict[str, Any] | None:
        normalized = {artifact_type.lower() for artifact_type in artifact_types}
        for artifact in artifacts:
            artifact_type = str(artifact.get("artifact_type", "")).lower()
            filename = str(artifact.get("filename", "")).lower()
            if artifact_type in normalized or filename in normalized:
                return artifact
        return None

    def _materialize_artifact(
        self,
        artifact: dict[str, Any],
        destination: Path,
        *,
        require_nonempty: bool = False,
        empty_error: str | None = None,
    ) -> None:
        destination.parent.mkdir(parents=True, exist_ok=True)

        for key in ("local_path", "source_path", "storage_key"):
            source_path = artifact.get(key)
            if not source_path:
                continue
            src = Path(str(source_path))
            if src.exists():
                shutil.copy(str(src), str(destination))
                self._ensure_materialized_file(destination, require_nonempty=require_nonempty, empty_error=empty_error)
                return

        download_url = artifact.get("download_url")
        if download_url:
            response = requests.get(str(download_url), timeout=60)
            response.raise_for_status()
            destination.write_bytes(response.content)
            self._ensure_materialized_file(destination, require_nonempty=require_nonempty, empty_error=empty_error)
            return

        raise FileNotFoundError(f"产物缺少可用源路径: {artifact}")

    def materialize_artifact(
        self,
        artifact: dict[str, Any],
        destination: Path,
        *,
        require_nonempty: bool = False,
        empty_error: str | None = None,
    ) -> None:
        self._materialize_artifact(
            artifact,
            destination,
            require_nonempty=require_nonempty,
            empty_error=empty_error,
        )

    @staticmethod
    def _ensure_materialized_file(
        destination: Path,
        *,
        require_nonempty: bool = False,
        empty_error: str | None = None,
    ) -> None:
        if not require_nonempty:
            return
        if destination.stat().st_size > 0:
            return
        raise ValueError(empty_error or f"{destination.name} 为空")
