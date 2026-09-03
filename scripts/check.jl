using Pkg

project_root = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(project_root)
Pkg.test()
