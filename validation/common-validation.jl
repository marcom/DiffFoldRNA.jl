using PythonCall

const EXIT_SUCCESS = 0
const EXIT_FAIL = 1

function parse_cmdline(args; n_start=2, n_end=2, atol=1e-7)
    print_usage() = println("usage: $PROGRAM_FILE [n_start] [n_end] [atol]")
    nargs = length(args)
    if any(a -> a in ("-h", "--help"), args)
        print_usage()
        return EXIT_FAIL, n_start, n_end, atol
    end
    if nargs == 1
        n_start = n_end = parse(Int, args[1])
    elseif nargs in (2, 3)
        n_start = parse(Int, args[1])
        n_end = parse(Int, args[2])
        if length(args) == 3
            atol = parse(Float64, args[3])
        end
    elseif nargs > 3
        println("error: too many arguments")
        print_usage()
        return EXIT_FAIL, n_start, n_end, atol
    end
    return EXIT_SUCCESS, n_start, n_end, atol
end

# should this go into julia module __init__() ?
function setup_python()
    pyloadpath = joinpath(@__DIR__, "..", "rna_prob_part_fn", "src")
    if ! isdir(pyloadpath)
        error("dir with python lib doesn't exist: $pyloadpath")
    end
    @py import sys
    sys.path.insert(0, pyloadpath)
    @info "py sys.path = $(sys.path)"
end
