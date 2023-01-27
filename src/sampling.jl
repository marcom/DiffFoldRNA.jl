using OffsetArrays

# Runs an O(n^3) precomp algorithm after which all structures are given an arbitrary unique ID.
# The structure for an ID can be retrieved in O(n^2) time.
# Useful for generating a uniform random sample from the space of structures for a sequence.
# Give a [None] sequence to allow any pairs.

# Uniform structure sampler
mutable struct AllStructs
    prim :: String
    dp :: OffsetMatrix{Int}
    hpmin :: Int
end
function AllStructs(prim::AbstractString, hpmin::Int)  # = AllStructs("", OffsetMatrix(zeros(0,0)))
    n = length(prim)
    # py: [[0]*len(prim) for _ in range(len(prim))]
    as = AllStructs(prim, OffsetMatrix(zeros(Int, n, n), 0:n-1, 0:n-1), hpmin)
    for i in axes(as.dp, 1)
        as.dp[i,i] = 1
    end
    HAIRPIN = hpmin
    for i in n-1:-1:0  # py: range(len(prim)-1, -1, -1):
        for j in i+1:n-1  # py: range(i+1, len(prim))
            as.dp[i, j] += as.dp[i+1, j]
            for k in i+HAIRPIN+1:j  # py: range(i+HAIRPIN+1, j+1)
                if ! valid_pair(as, i, k)
                    continue
                end
                # TODO: ensure branchlessness
                as.dp[i, j] += (i+1 < k-1 ? as.dp[i+1, k-1] : 1) * (k+1 < j ? as.dp[k+1, j] : 1)
            end
        end
    end
    return as
end
AllStructs(n::Int, hpmin::Int) = AllStructs("N"^n, hpmin)

function valid_pair(as::AllStructs, i::Int, j::Int)
    # TODO: this is because the dp-matrix as.dp indexes from 0, and
    # String indexes from 1
    i += 1
    j += 1

    # Always allow 'N' to pair with anything.
    if as.prim[i] == 'N' || as.prim[j] == 'N'
        return true
    end
    return valid_pair(as.prim[i], as.prim[j])
end

count_structures(as::AllStructs) = as.dp[begin, end]

# Gets the nth structure where n is in [0, self.count_structures()]
function get_nth(as::AllStructs, n::Int)
    len = length(as.prim)
    match = [i for i in 1:len]  # OffsetVector([i for i in 1:len], 0:len-1)
    HAIRPIN = as.hpmin

    function trace(i, j, n)
        if i >= j
            return
        end
        if n < as.dp[i+1, j]
            trace(i+1, j, n)
            return
        end
        n -= as.dp[i+1, j]
        for k in i+HAIRPIN+1:j  # py: range(i+HAIRPIN+1, j+1):
            if ! valid_pair(as, i, k)
                continue
            end
            left = i+1 < k-1 ? as.dp[i+1, k-1] : 1
            right = k+1 < j ? as.dp[k+1, j] : 1
            if n < left * right
                # TODO: indices adjusted because match is 1-based
                # old
                # match[i] = k
                # match[k] = i
                match[i+1] = k+1
                match[k+1] = i+1
                trace(i+1, k-1, fld(n, right))  # py: n//right)
                trace(k+1, j, rem(n, right))  # py: n % right)
                return
            end
            n -= left * right
        end
    end

    trace(0, len-1, n)
    return match
end

random_matching(as::AllStructs) = get_nth(as, rand(0:count_structures(as)))
random_dbn(as::AllStructs) = matching_to_dbn(random_matching(as))
