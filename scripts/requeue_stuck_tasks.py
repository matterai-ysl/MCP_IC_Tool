from src.vasp_server.task_manager.manager import TaskManager


def requeue_stuck_tasks(task_manager: TaskManager | None = None) -> dict[str, int]:
    manager = task_manager or TaskManager()
    expired_leases = manager.requeue_expired_leases()
    orphaned_running = manager.mark_orphaned_running_tasks()
    return {
        "expired_leases": expired_leases,
        "orphaned_running": orphaned_running,
    }


def main() -> int:
    result = requeue_stuck_tasks()
    print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
