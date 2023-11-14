using Plots
gr()
data = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0;
        0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0;
        0 0 0 0 0 1 1 0.5 0.5 1 1 0 0 0 0 0;
        0 0 0 0 0 1 1 0.5 0.5 1 1 0 0 0 0 0;
        0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0;
        0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;]
hm = heatmap(1:size(data,1),
    1:size(data,2), data,
    # c=cgrad([:white,:red]),
    c = cgrad([:white,:green,:red], [0, 1, 0.5], categorical = false), 
    xticks = 0:1:size(data,1),
    yticks = 0:1:size(data,2),
    aspect_ratio=1, 
    size=(600, 600),
    yflip = true,
    xmirror = true,
    colorbar = false,
    xlabel = "y-axis",
    ylabel = "x-axis",
    title = "z = 9"
    )

    m = size(data,1) 
    n = size(data,2)
    vline!(0.5:(n+0.5), c=:black, legend=:none)
    hline!(0.5:(m+0.5), c=:black, legend=:none)
    
    ann = [(i,j,0) for i in 1:5 for j in 1:16]
    annotate!(ann)
    ann = [(i,j,0) for i in 6:11 for j in 1:5]
    annotate!(ann)
    ann = [(i,j,0) for i in 6:11 for j in 12:16]
    annotate!(ann)
    ann = [(i,j,0) for i in 12:16 for j in 1:16]
    annotate!(ann)
    annotate!((7,7,1))
    annotate!((7,10,1))
    annotate!((10,10,1))
    annotate!((10,7,1))
    annotate!((7,8,1))
    annotate!((7,9,1))
    annotate!((8,7,1))
    annotate!((8,8,-15))
    annotate!((8,9,-15))
    annotate!((9,7,1))
    annotate!((9,8,-15))
    annotate!((9,9,-15))
    annotate!((8,10,1))
    annotate!((9,10,1))
    annotate!((10,8,1))
    annotate!((10,9,1))
    annotate!((6,6,0))
    annotate!((6,7,1))
    annotate!((7,6,1))
    annotate!((6,10,1))
    annotate!((6,11,0))
    annotate!((7,11,1))
    annotate!((10,6,1))
    annotate!((11,6,0))
    annotate!((11,7,1))
    annotate!((10,11,1))
    annotate!((11,11,0))
    annotate!((11,10,1))
    annotate!((6,8,1))
    annotate!((6,9,1))
    annotate!((8,6,1))
    annotate!((9,6,1))
    annotate!((11,8,1))
    annotate!((11,9,1))
    annotate!((9,11,1))
    annotate!((8,11,1))


