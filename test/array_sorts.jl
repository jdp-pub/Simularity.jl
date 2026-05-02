"""
    bubble_sort(A::AbstractArray,k::AbstractFloat=1,mode::String="ascending")

# Arguments
- `A::AbstractArray`: An array of objects (needs to be tested for lexicographical support)
- `k::Float`: The percentage of greatest elements to gurantee gets sorted.
- `mode::String` : "AS" ascending or "DS" descending. All other arguments default to ascending.

# Return 
A fully sorted or partially sorted array.

# Description
Bubble sort is a simple array sorting algorithm[^bubble-sort] with the property 
that larger values are guaranteed to arive at the and of the list 
each pass (hence the bubbling). This leads to (some interesting implications[^when-would-you-ever-want-bubblesort] 
and is best used as a way to quickly partially order lists, most relevant in 
certain simulation tasks[^depth-sorting-of-billboard-particles-how-can-i-do-it]

# References
[^bubble-sort]: [Bubble Sort, https://en.wikipedia.org/wiki/Bubble_sort (accessed 4 14, 2026).](https://en.wikipedia.org/wiki/Bubble_sort)
[^when-would-you-ever-want-bubblesort]: [When Would You Ever Want Bubble Sort?, https://buttondown.com/hillelwayne/archive/when-would-you-ever-want-bubblesort/ (accessed 4 14, 2026).](https://buttondown.com/hillelwayne/archive/when-would-you-ever-want-bubblesort/)
[^depth-sorting-of-billboard-particles-how-can-i-do-it]: [Depth Sorting of Billboard Particles, How Can I Do It?, https://discussions.unity.com/t/depth-sorting-of-billboard-particles-how-can-i-do-it/5053 (accessed 4 14, 2026).](https://discussions.unity.com/t/depth-sorting-of-billboard-particles-how-can-i-do-it/5053). 
"""
function bubble_sort(A::AbstractArray;k::AbstractFloat=1.,mode::String="AS")
    V = A

    sorted = false 
    swaps = 0
    vtemp = 0
    vtemp = typeof(A[1])

    kn = Int(length(A)*k)
    kx = 0

    (mode == "DS") ? op = (<) : op = (>)

    while !sorted
        swaps = 0
        for nx in 1:length(V)-1
            if op(V[nx],V[nx+1])
                vtemp = V[nx]
                V[nx] = V[nx+1]
                V[nx+1] = vtemp
                swaps = swaps+1
            end
        end

        kx = kx + 1
        if swaps == 0 || kx == kn
            sorted = true
        end
    end
    return V 
end

"""
    insertion_sort(A::AbstractArray)

# Arguments
- `A::AbstractArray`: An array of objects (needs to be tested for lexicographical support)
- `mode::String` : "AS" ascending or "DS" descending. All other arguments default to ascending.

# Return 
A sorted array.

# Description

# References
"""
function insertion_sort(A::AbstractArray;mode::String="AS")
    V = A
    nt = 0
    xtemp = 0

    (mode == "DS") ? op = (<) : op = (>)

    for nx = 1:length(V)-1
        if op(V[nx],V[nx+1])
            nt = nx
            while nt >=1 && op(V[nt],V[nt+1])
                xtemp = V[nt]
                V[nt] = V[nt+1]
                V[nt+1] = xtemp
                nt = nt-1
            end
        end
    end
    return V
end

"""
    merge_sort(A::AbstractArray)

# Arguments
- `A::AbstractArray`: An array of objects (needs to be tested for lexicographical support).
- `mode::String` : "AS" ascending or "DS" descending. All other arguments default to ascending.

# Return 
A sorted array.

# Description

# References
"""
function merge_sort(A::AbstractArray;mode::String="AS")
    if length(A) == 1
        return A
    end

    V1 = merge_sort(A[1:Int(floor(length(A)/2))],mode=mode)
    V2 = merge_sort(A[Int(floor(length(A)/2))+1:length(A)],mode=mode)

    V = vcat(V1,V2)
    return insertion_sort(V,mode=mode)
end
