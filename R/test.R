fun <- function(t) {
    ap <- lapply(t, unique)
    cd <- ap == 1
    ll <- lengths(cd)
    all(ll)
}
