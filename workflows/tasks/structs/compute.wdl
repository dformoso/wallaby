version 1.0

struct Resources {
    Int cpu
    String memory_gb
    String disks
    String zones
    Int preemptible
    Int maxRetries
}

struct Compute {
    Map[String, Resources] size
}