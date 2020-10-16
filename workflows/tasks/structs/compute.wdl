version 1.0

struct Resources {
    Int cpu
    Int memory_gb
    String disks
    String zones
    Int preemptible
    Int maxRetries
    String gpuType
    Int gpuCount
}

struct Compute {
    Map[String, Resources] size
}