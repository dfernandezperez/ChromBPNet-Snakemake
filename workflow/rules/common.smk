# --- Partition of the GPU jobs ---
# irb_gpu_3090 runs at most `gpu_lane_limit` GPUs per user at a time; spot_gpu has no such limit, but its
# jobs can be preempted (they are requeued and restart from scratch, see `slurm-requeue` in the profile).
# A job can't be submitted to both partitions (they accept different QOS), so the partition of each GPU
# job is chosen when snakemake submits it, from the rule's `slurm_partition` in resources.yaml:
#   main: always the main partition (model training: long jobs that should not be preempted)
#   auto: the main partition while the user has fewer than `gpu_lane_limit` jobs there (running or
#         pending), otherwise the overflow partition. Retries of failed jobs always go to the main one.
#   any other value: that partition
# Decisions are logged to logs/gpu_routing.log.
import getpass
import subprocess
import threading
import time

GPU_ROUTING = RESOURCES["gpu_routing"]

# Jobs that become ready together are all evaluated before any of them is submitted, so the decisions
# taken since the last squeue snapshot are counted too. The snapshot is only refreshed between
# scheduling rounds, once the previous decisions have been submitted and appear in squeue.
_lanes = {"lock": threading.Lock(), "snapshot_time": float("-inf"), "snapshot_jobs": 0,
          "claimed_since_snapshot": 0, "last_claim": float("-inf"), "decisions": {}}


def _main_partition_jobs():
    """Number of the user's jobs (pending or running) in the main GPU partition."""
    out = subprocess.run(
        ["squeue", "-h", "-u", getpass.getuser(), "-p", GPU_ROUTING["main_partition"],
         "-t", "PENDING,RUNNING,CONFIGURING,COMPLETING,SUSPENDED,REQUEUED,REQUEUE_HOLD", "-o", "%i"],
        capture_output=True, text=True, timeout=30, check=True).stdout
    return len(out.split())


def _decision_key(rulename, wildcards, attempt):
    return (rulename, tuple(sorted(dict(wildcards).items())), attempt)


def _waiting_main_jobs():
    """Ready jobs of `main` rules (model training) that have not been routed yet.

    The jobs that become ready together are routed in an arbitrary order, so the main-partition lanes
    they will need are reserved for them before routing `auto` jobs there.
    """
    try:
        return sum(1 for job in workflow.dag.ready_jobs
                   if RESOURCES.get(job.rule.name, {}).get("slurm_partition") == "main"
                   and _decision_key(job.rule.name, job.wildcards_dict, job.attempt) not in _lanes["decisions"])
    except Exception:
        return 0


def _choose_partition(wildcards, attempt, rulename, always_main):
    main, overflow = GPU_ROUTING["main_partition"], GPU_ROUTING["overflow_partition"]
    # The Snakefile is parsed again inside every cluster job, where the partition is irrelevant
    if workflow.dryrun or not workflow.is_main_process:
        return main
    always_main = always_main or not GPU_ROUTING["enabled"] or attempt > 1
    key = _decision_key(rulename, wildcards, attempt)
    now = time.time()
    with _lanes["lock"]:
        # In case the resources of a job are evaluated more than once
        decision = _lanes["decisions"].get(key)
        if decision is not None and now - decision[1] < GPU_ROUTING["claim_ttl_seconds"]:
            return decision[0]
        if (now - _lanes["snapshot_time"] > GPU_ROUTING["squeue_ttl_seconds"]
                and now - _lanes["last_claim"] > 2):
            try:
                _lanes["snapshot_jobs"] = _main_partition_jobs()
            except Exception:
                _lanes["snapshot_jobs"] = GPU_ROUTING["gpu_lane_limit"] # Unknown: assume no free lane
            _lanes["snapshot_time"], _lanes["claimed_since_snapshot"] = now, 0
        used = _lanes["snapshot_jobs"] + _lanes["claimed_since_snapshot"]
        if not always_main:
            used += _waiting_main_jobs()
        partition = main if always_main or used < GPU_ROUTING["gpu_lane_limit"] else overflow
        if partition == main:
            _lanes["claimed_since_snapshot"] += 1
        _lanes["decisions"][key] = (partition, now)
        _lanes["last_claim"] = now
    os.makedirs("logs", exist_ok=True)
    with open("logs/gpu_routing.log", "a") as log:
        log.write(f"{time.strftime('%Y-%m-%d %H:%M:%S')}\t{rulename}\t{dict(wildcards)}\t"
                  f"attempt={attempt}\tmain_partition_jobs={used}\t{partition}\n")
    return partition


# The `input` argument is unused, but it makes snakemake skip these functions when it first evaluates
# the resources of every job (while building the DAG) and call them only once the job is ready to run.
def _route_auto(wildcards, input, attempt, rulename):
    return _choose_partition(wildcards, attempt, rulename, always_main=False)


def _route_main(wildcards, input, attempt, rulename):
    # Counted as a claim of a main-partition lane, so that `auto` jobs evaluated in the same
    # scheduling round don't take it
    return _choose_partition(wildcards, attempt, rulename, always_main=True)


def gpu_partition(rule_name):
    """`slurm_partition` resource of a GPU rule (see the explanation above)."""
    partition = RESOURCES[rule_name]["slurm_partition"]
    return {"auto": _route_auto, "main": _route_main}.get(partition, partition)


# --- Thread limits ---
# TensorFlow, numpy (OpenBLAS) and numba size their thread pools with the cores of the whole node
# instead of the CPUs allocated to the job. Prefix commands with this to use `n` threads.
def thread_env(n):
    n = max(1, int(n))
    return (f"OMP_NUM_THREADS={n} OPENBLAS_NUM_THREADS={n} MKL_NUM_THREADS={n} NUMBA_NUM_THREADS={n} "
            f"TF_NUM_INTRAOP_THREADS={n} TF_NUM_INTEROP_THREADS=1")
