module Ajna

#   This file is part of Ajna.jl. It is licensed under the AGPL license
#   Ajna Copyright (C) 2020 Michael Reed

using Grassmann, ColorTypes, ImageInTerminal

struct Rectangle
    xmin::Float64
    xmax::Float64
    ymin::Float64
    ymax::Float64
    nx::Int
    ny::Int
end

Rectangle(x) = Rectangle(-x,x,-x,x,100,100)

function raster(ga::Vector,R::Rectangle=Rectangle(3))
    out = zeros(GrayA{Float64},R.nx,R.ny)
    δx,δy = (R.xmax-R.xmin)/(R.nx-1),(R.ymax-R.ymin)/(R.ny-1)
    δ = sqrt(δx^2+δy^2)/2
    for x ∈ 1:R.nx, y ∈ 1:R.ny
        P = Chain(1.0,R.xmin+(x-1)*δx,R.ymin+(y-1)*δy)
        c = 0.0
        for g ∈ ga
            norm(P∧g)<δ && (c+=1)
        end
        out[x,y] = GrayA(c,c)
    end
    return out
end

end # module
