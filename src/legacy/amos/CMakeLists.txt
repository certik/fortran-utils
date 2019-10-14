set(SRC
    d1mach.f90
    i1mach.f90
    r1mach.f90
    dgamln.f
    dsclmr.f
    fdump.f
    zabs.f
    zacai.f
    zacon.f
    zairy.f
    zasyi.f
    zbesh.f
    zbesi.f
    zbesj.f
    zbesk.f
    zbesy.f
    zbinu.f
    zbiry.f
    zbknu.f
    zbuni.f
    zbunk.f
    zdiv.f
    zexp.f
    zkscl.f
    zlog.f
    zmlri.f
    zmlt.f
    zrati.f
    zs1s2.f
    zseri.f
    zshch.f
    zsqrt.f
    zuchk.f
    zunhj.f
    zuni1.f
    zuni2.f
    zunik.f
    zunk1.f
    zunk2.f
    zuoik.f
    zwrsk.f
)

add_library(amos ${SRC})
target_compile_options(amos PRIVATE -Wno-implicit-interface -Wno-compare-reals)
add_library(p::amos ALIAS amos)
