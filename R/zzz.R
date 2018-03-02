## This file is absolutely essential. Without it, R won't load your module and won't know what you're talking about
loadModule("RmmquantParameters_module", TRUE)

.onUnload <- function (libpath) {
    library.dynam.unload("Rmmquant", libpath)
}
