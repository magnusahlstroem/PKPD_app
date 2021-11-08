get_named_args <- function () {
    as.list( match.call(
        def = sys.function( -1 ),
        call = sys.call(-1)) )[-1]
}
