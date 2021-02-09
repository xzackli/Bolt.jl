using BinDeps

BinDeps.@setup

librecfast = library_dependency("librecfast")

provides(SimpleBuild,
    (@build_steps begin
        ChangeDirectory(joinpath(BinDeps.depsdir(librecfast),"src","recfast"))
        `mkdir -p $(libdir(librecfast))`
        `gfortran -ffixed-line-length-none -fPIC -shared -g -O0 recfast.for -o $(libdir(librecfast))/librecfast.so`
    end), librecfast, os = :Unix)
    
BinDeps.@install Dict(:librecfast=>:librecfast)
