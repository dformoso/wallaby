##################################
########      STRUCTS     ########
##################################
version development

struct BWAIndex {
    File fasta
    File amb
    File ann
    File bwt
    File pac
    File sa
}

struct SplitBAMs {
    File MM
    File MU
    File UM
    File UU
}

struct SplitQNAMEs {
    File MM
    File MU
    File UM
    File UU
}