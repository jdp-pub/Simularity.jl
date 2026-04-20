using BenchmarkTools

"""
    f_compare(func_list::AbstractArray,n::AbstractArray,avg_n::Int=10)


# Arguments
- `func_list::AbstractArray`: An array of functions to compare. 
- `n::AbstractArray`: An array of arguments for the functions to compare.
- `avg_n::Int`: The number of runs per argument, suppresses numerical noise.

# Return 
An array of times corresponding to each function.

# Description
For benchmarking similar algorithm run times.

"""
function f_compare(func_list::AbstractArray,n::AbstractArray,avg_n::Int=10)

    times = zeros(length(func_list))
    time_avg = zeros(Float64,avg_n)

    for fx in 1:length(func_list)
        for tx in 1:avg_n
            time_avg[tx] = @belapsed $func_list[fx]($n[fx])
        end
        times[fx] = mean(time_avg)
    end
    return times
end

