# ---------------------------------------------------------------------------- #
Parse_Arguments = function(
    name = NULL
){
    if(is.null(name))
        stop('Please use a valid function name')

    library(argparser, quietly=TRUE)

    p = arg_parser(name)
    p = add_argument(p, "--basedir", help="basedir", type="list", default = NULL)
    p = add_argument(p, "--input",   help="input",   type="list", default = NULL)
    p = add_argument(p, "--output",  help="output",  type="list", default = NULL)
    p = add_argument(p, "--params",  help="params",  type="list", default = NULL)

    # argv = parse_args(p)
    argv = parse_args(p)[]
    snames = c('input','output','params')
    argv = lapply(setNames(snames, snames), function(x)jsonlite::fromJSON(argv[[x]][[1]]))
    return(argv)
}
