# Plate geometry.
#
# SBS-format plates share one footprint regardless of density: the well field is
# 72 x 108 mm whether the plate is 96-well or 384-well. Only the pitch changes.
# Fitting the GP in millimetres therefore makes a lengthscale prior portable
# across plate formats, where a prior in well indices is not.

const SBS_FIELD_MM = (rows = 72.0, cols = 108.0)

"""
    well_pitch(nrow, ncol) -> Float64

Centre-to-centre well spacing in mm for an SBS plate of the given shape:
9.0 mm for 8x12, 4.5 mm for 16x24, 2.25 mm for 32x48.
"""
function well_pitch(nrow::Integer, ncol::Integer)
    p = SBS_FIELD_MM.rows / nrow
    q = SBS_FIELD_MM.cols / ncol
    isapprox(p, q; rtol = 1e-6) ||
        error("well_pitch: $(nrow)x$(ncol) is not an SBS layout (row pitch $p mm, column pitch $q mm)")
    return p
end
