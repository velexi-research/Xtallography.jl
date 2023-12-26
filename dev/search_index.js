var documenterSearchIndex = {"docs":
[{"location":"functions/","page":"Functions","title":"Functions","text":"CurrentModule = XtallographyUtils","category":"page"},{"location":"functions/#Functions","page":"Functions","title":"Functions","text":"","category":"section"},{"location":"functions/","page":"Functions","title":"Functions","text":"","category":"page"},{"location":"functions/#Lattice-System-Functions","page":"Functions","title":"Lattice System Functions","text":"","category":"section"},{"location":"functions/","page":"Functions","title":"Functions","text":"is_bravais_lattice","category":"page"},{"location":"functions/#XtallographyUtils.is_bravais_lattice","page":"Functions","title":"XtallographyUtils.is_bravais_lattice","text":"is_bravais_lattice(lattice_system::LatticeSystem, centering::Centering) -> Bool\n\nis_bravais_lattice(unit_cell::UnitCell) -> Bool\n\nDetermine if the unit cell defined by unit_cell or lattice_system and centering is a valid Bravais lattice type.\n\nReturn values\n\ntrue if lattice_system and centering define a valid Bravais lattice type; false otherwise\n\nExamples\n\njulia> is_bravais_lattice(cubic, body_centered)\ntrue\n\njulia> is_bravais_lattice(cubic, base_centered)\nfalse\n\njulia> is_bravais_lattice(UnitCell(TetragonalLatticeConstants(2, 3), primitive))\ntrue\n\njulia> is_bravais_lattice(UnitCell(TetragonalLatticeConstants(2, 3), face_centered))\nfalse\n\n\n\n\n\n","category":"function"},{"location":"functions/","page":"Functions","title":"Functions","text":"","category":"page"},{"location":"functions/#Unit-Cell-Functions","page":"Functions","title":"Unit Cell Functions","text":"","category":"section"},{"location":"functions/#Unit-Cell-Properties","page":"Functions","title":"Unit Cell Properties","text":"","category":"section"},{"location":"functions/","page":"Functions","title":"Functions","text":"basis\nsurface_area(::UnitCell)\nvolume(::UnitCell)","category":"page"},{"location":"functions/#XtallographyUtils.basis","page":"Functions","title":"XtallographyUtils.basis","text":"basis(unit_cell::UnitCell) -> (Vector{Float64}, Vector{Float64}, Vector{Float64})\n\nbasis(\n    lattice_constants::LatticeConstants\n) -> (Vector{Float64}, Vector{Float64}, Vector{Float64})\n\nReturn a set of basis vectors veca vecb vecc for the unit cell defined by unit_cell or lattice_constants.\n\nReturn values\n\nbasis vectors veca, vecb, vecc\n\nExamples\n\njulia> B = basis(UnitCell(LatticeConstants([1, 0, 0], [1, 1, 0], [1, 0, 2]), primitive));\n\njulia> B[1] ≈ [1.0, 0.0, 0.0]\ntrue\njulia> B[2] ≈ [1.0, 1.0, 0.0]\ntrue\njulia> B[3] ≈ [1.0, 0.0, 2.0]\ntrue\n\njulia> B = basis(LatticeConstants([1, 0, 0], [1, 1, 0], [1, 0, 2]));\n\njulia> B[1] ≈ [1.0, 0.0, 0.0]\ntrue\njulia> B[2] ≈ [1.0, 1.0, 0.0]\ntrue\njulia> B[3] ≈ [1.0, 0.0, 2.0]\ntrue\n\n\n\n\n\n","category":"function"},{"location":"functions/#XtallographyUtils.surface_area-Tuple{UnitCell}","page":"Functions","title":"XtallographyUtils.surface_area","text":"surface_area(unit_cell::UnitCell) -> Float64\n\nsurface_area(lattice_constants::LatticeConstants) -> Float64\n\nCompute the surface area of the unit cell defined by unit_cell or lattice_constants.\n\nReturn values\n\nsurface area of the unit cell\n\nExamples\n\njulia> surface_area(UnitCell(LatticeConstants([1, 0, 0], [1, 1, 0], [1, 0, 1]), primitive))\n7.464101615137754\n\njulia> surface_area(LatticeConstants([1, 0, 0], [1, 1, 0], [1, 0, 1]))\n7.464101615137754\n\n\n\n\n\n","category":"method"},{"location":"functions/#XtallographyUtils.volume-Tuple{UnitCell}","page":"Functions","title":"XtallographyUtils.volume","text":"volume(unit_cell::UnitCell) -> Float64\n\nvolume(lattice_constants::LatticeConstants) -> Float64\n\nCompute the volume of the unit cell defined by unit_cell or lattice_constants.\n\nReturn values\n\nvolume of the unit cell\n\nExamples\n\njulia> volume(UnitCell(LatticeConstants([1, 0, 0], [1, 1, 0], [1, 0, 2]), primitive))\n2.0\n\njulia> volume(LatticeConstants([1, 0, 0], [1, 1, 0], [1, 0, 2]))\n2.0\n\n\n\n\n\n","category":"method"},{"location":"functions/#Unit-Cell-Standardization-and-Comparison-Functions","page":"Functions","title":"Unit Cell Standardization and Comparison Functions","text":"","category":"section"},{"location":"functions/","page":"Functions","title":"Functions","text":"conventional_cell\nis_equivalent_unit_cell\nis_supercell\nreduced_cell\nstandardize(::UnitCell)\nstandardize(::LatticeConstants, ::Centering)\nstandardize(::LatticeConstants)","category":"page"},{"location":"functions/#XtallographyUtils.conventional_cell","page":"Functions","title":"XtallographyUtils.conventional_cell","text":"conventional_cell(unit_cell::UnitCell) -> UnitCell\n\nReturn the IUCr conventional cell that is equivalent to unit_cell.\n\nReturn values\n\nIUCr conventional cell for unit_cell\n\nExamples\n\njulia> unit_cell = UnitCell(LatticeConstants([1, 1, 0], [1, -1, 0], [0, 1, 1]), primitive);\n\njulia> lattice_system(unit_cell.lattice_constants)\nTriclinic()\n\njulia> conventional_cell(unit_cell) ≈ UnitCell(CubicLatticeConstants(2.0), FaceCentered())\ntrue\n\n\n\n\n\n","category":"function"},{"location":"functions/#XtallographyUtils.is_equivalent_unit_cell","page":"Functions","title":"XtallographyUtils.is_equivalent_unit_cell","text":"is_equivalent_unit_cell(\n    unit_cell_test::UnitCell,\n    unit_cell_ref::UnitCell;\n    tol::Real=1e-3\n) -> Bool\n\nCheck if the unit cell defined by unit_cell_test is equivalent to the unit cell defined by unit_cell_ref.\n\nKeyword Arguments\n\ntol: absolute tolerance of the deviation between the reduced unit cells defined by unit_cell_test and unit_cell_ref.\n\nReturn values\n\ntrue if the test unit cell is equivalent to the reference unit cell; false otherwise\n\nExamples\n\njulia> lattice_constants_ref = LatticeConstants([1, 0, 0], [1, 1, 0], [0, 0, 2]);\n\njulia> unit_cell_ref = UnitCell(lattice_constants_ref, body_centered);\n\njulia> lattice_constants_test = TetragonalLatticeConstants(1.0, 2.0);\n\njulia> unit_cell_test = UnitCell(lattice_constants_test, body_centered);\n\njulia> is_equivalent_unit_cell(unit_cell_test, unit_cell_ref)\ntrue\n\n\n\n\n\nis_equivalent_unit_cell(\n    lattice_constants_test::LatticeConstants,\n    lattice_constants_ref::LatticeConstants;\n    tol::Real=1e-3\n) -> Bool\n\nCheck if the primitive unit cell defined by lattice_constants_test is equivalent to the unit cell defined by lattice_constants_ref.\n\nKeyword Arguments\n\ntol: absolute tolerance of the deviation between the reduced unit cells defined by lattice_constants_test and lattice_constants_ref.\n\nReturn values\n\ntrue if the test unit cell is equivalent to the reference unit cell; false otherwise\n\nExamples\n\njulia> lattice_constants_ref = LatticeConstants([1, 0, 0], [1, 1, 0], [0, 0, 2]);\n\njulia> is_equivalent_unit_cell(TetragonalLatticeConstants(1.0, 2.0), lattice_constants_ref)\ntrue\n\n\n\n\n\n","category":"function"},{"location":"functions/#XtallographyUtils.is_supercell","page":"Functions","title":"XtallographyUtils.is_supercell","text":"is_supercell(\n    lattice_constants_test::LatticeConstants,\n    lattice_constants_ref::LatticeConstants;\n    tol::Real=1e-3,\n    max_index::Integer=3\n) -> Bool\n\nCheck if the unit cell defined by lattice_constants_test is a supercell of the unit cell defined by lattice_constants_ref.\n\nKeyword Arguments\n\ntol: absolute tolerance of the deviation between test unit cell and supercells of the reference unit cell\nmax_index: maximum multiple of the basis vectors of the reference unit cell to include when checking whether the test unit cell is a supercell of the reference unit cell\nnote: Note\nmax_index is ignored for Cubic lattice systems.\n\nReturn values\n\ntrue if the test unit cell is a supercell of the reference unit cell; false otherwise\n\nExamples\n\njulia> lattice_constants_ref = LatticeConstants([1, 0, 0], [0, 1, 0], [0, 0, 1]);\n\njulia> is_supercell(CubicLatticeConstants(2), lattice_constants_ref)\ntrue\n\njulia> is_supercell(CubicLatticeConstants(2.5), lattice_constants_ref)\nfalse\n\njulia> is_supercell(LatticeConstants([1, 0, 0], [0, 2, 0], [0, 0, 3]),\n                    lattice_constants_ref)\nfalse\n\n\n\n\n\n","category":"function"},{"location":"functions/#XtallographyUtils.reduced_cell","page":"Functions","title":"XtallographyUtils.reduced_cell","text":"reduced_cell(unit_cell::UnitCell) -> UnitCell\n\nCompute the primitive reduced cell for the lattice defined by unit_cell. The Selling-Delaunay reduction algorithm is used to compute the reduced basis.\n\nReturn values\n\nprimitive reduced cell\n\nExamples\n\njulia> reduced_cell(UnitCell(LatticeConstants([1, 0, 0], [1, 1, 0], [0, 0, 2]), primitive))\nUnitCell(TetragonalLatticeConstants(1.0, 2.0), Primitive())\n\n\n\n\n\n","category":"function"},{"location":"functions/#XtallographyUtils.standardize-Tuple{UnitCell}","page":"Functions","title":"XtallographyUtils.standardize","text":"standardize(unit_cell::UnitCell) -> LatticeConstants\n\nStandardize the lattice constants and centering for unit_cell.\n\nnote: Note\nThis function only enforces lattice constant constraints. It does not modify the Bravais lattice type. To find an equivalent Bravais lattice with higher symmetry (if one exists), use conventional_cell().\n\nnote: Note\nLattice constant standardizations are based on the conventions provided in the Table 3.1.4.1. of the International Tables for Crystallography (2016). For triclinic lattices, the lattice constants are standardized using the following conventions:a ≤ b ≤ c\nall three angles are acute (Type I cell) or all three angles are non-acute (Type II cell)\nangles sorted in increasing order when edge lengths are equal\nα ≤ β when a = b\nβ ≤ γ when b = c\nα ≤ β ≤ γ when a = b = c\n\nnote: Note\nExcept for monoclinic lattices, standardize() does not modify the unit cell.For monoclinic lattices, the unit cell may be potentially modified in two ways:the 2D unit cell in the plane normal to the b-axis may be reduced in order to satisfy the IUCr conventions for a, c, and β;\nbase-centered unit cells are transformed to equivalent body-centered unit cells.\n\nReturn values\n\nunit cell with standardized lattice constants and centering\n\nExamples\n\njulia> unit_cell = UnitCell(OrthorhombicLatticeConstants(3, 2, 1), primitive)\nUnitCell(OrthorhombicLatticeConstants(3.0, 2.0, 1.0), Primitive())\njulia> standardize(unit_cell)\nUnitCell(OrthorhombicLatticeConstants(1.0, 2.0, 3.0), Primitive())\n\njulia> unit_cell = UnitCell(OrthorhombicLatticeConstants(3, 2, 1), base_centered)\nUnitCell(OrthorhombicLatticeConstants(3.0, 2.0, 1.0), BaseCentered())\njulia> standardize(unit_cell)\nUnitCell(OrthorhombicLatticeConstants(2.0, 3.0, 1.0), BaseCentered())\n\n\n\n\n\n","category":"method"},{"location":"functions/#XtallographyUtils.standardize-Tuple{LatticeConstants, Centering}","page":"Functions","title":"XtallographyUtils.standardize","text":"standardize(\n    lattice_constants::LatticeConstants, centering::Centering\n) -> (LatticeConstants, Centering)\n\nStandardize the lattice constants and centering for the unit cell defined by lattice_constants and centering.\n\nReturn values\n\nstandardized lattice constants and centering\n\nExamples\n\njulia> standardize(OrthorhombicLatticeConstants(3, 2, 1), primitive)\n(OrthorhombicLatticeConstants(1.0, 2.0, 3.0), Primitive())\n\njulia> standardize(OrthorhombicLatticeConstants(3, 2, 1), base_centered)\n(OrthorhombicLatticeConstants(2.0, 3.0, 1.0), BaseCentered())\n\n\n\n\n\n","category":"method"},{"location":"functions/#XtallographyUtils.standardize-Tuple{LatticeConstants}","page":"Functions","title":"XtallographyUtils.standardize","text":"standardize(lattice_constants::LatticeConstants) -> LatticeConstants\n\nStandardize the lattice constants the primitive unit cell defined by lattice_constants.\n\nReturn values\n\nstandardized lattice constants for primitive unit cell\n\nExamples\n\njulia> standardize(OrthorhombicLatticeConstants(3, 2, 1))\nOrthorhombicLatticeConstants(1.0, 2.0, 3.0)\n\n\n\n\n\n","category":"method"},{"location":"functions/#Other-Functions","page":"Functions","title":"Other Functions","text":"","category":"section"},{"location":"functions/","page":"Functions","title":"Functions","text":"lattice_system","category":"page"},{"location":"functions/#XtallographyUtils.lattice_system","page":"Functions","title":"XtallographyUtils.lattice_system","text":"lattice_system(lattice_constants::LatticesConstants) -> LatticeSystem\n\nReturn the lattice system for a set of lattice_constants.\n\nReturn values\n\nlattice system\n\nExamples\n\njulia> lattice_system(HexagonalLatticeConstants(2, 4))\nHexagonal()\n\njulia> lattice_system(CubicLatticeConstants(2))\nCubic()\n\n\n\n\n\n","category":"function"},{"location":"functions/#Lattice-Specific-Functions","page":"Functions","title":"Lattice-Specific Functions","text":"","category":"section"},{"location":"functions/#Triclinic-Systems","page":"Functions","title":"Triclinic Systems","text":"","category":"section"},{"location":"functions/","page":"Functions","title":"Functions","text":"convert_to_mP\nconvert_to_mI\nconvert_to_mC\nis_triclinic_type_I_cell\nsatisfies_triclinic_angle_constraints","category":"page"},{"location":"functions/#XtallographyUtils.convert_to_mP","page":"Functions","title":"XtallographyUtils.convert_to_mP","text":"convert_to_mP(\n    lattice_constants::TriclinicLatticeConstants\n) -> MonoclinicLatticeConstants\n\nAttempt to convert the triclinic unit cell defined by lattice_constants to an equivalent primitive monoclinic unit cell.\n\nReturn values\n\nlattice constants for the equivalent primitive monoclinic unit cell if one exists\n\nExceptions\n\nThrows an ErrorException if the triclinic unit cell defined by lattice_constants is not equivalent to a primitive monoclinic unit cell.\n\n\n\n\n\n","category":"function"},{"location":"functions/#XtallographyUtils.convert_to_mI","page":"Functions","title":"XtallographyUtils.convert_to_mI","text":"convert_to_mI(\n    lattice_constants::TriclinicLatticeConstants\n) -> MonoclinicLatticeConstants\n\nAttempt to convert the triclinic unit cell defined by lattice_constants to an equivalent body-centered monoclinic unit cell.\n\nReturn values\n\nlattice constants for the equivalent body-centered monoclinic unit cell if one exists; nothing otherwise\n\nExceptions\n\nThrows an ErrorException if the triclinic unit cell defined by lattice_constants is not equivalent to a body-centered monoclinic unit cell.\n\n\n\n\n\n","category":"function"},{"location":"functions/#XtallographyUtils.convert_to_mC","page":"Functions","title":"XtallographyUtils.convert_to_mC","text":"convert_to_mC(\n    lattice_constants::TriclinicLatticeConstants\n) -> MonoclinicLatticeConstants\n\nAttempt to convert the triclinic unit cell defined by lattice_constants to an equivalent base-centered monoclinic unit cell.\n\nReturn values\n\nlattice constants for the equivalent base-centered monoclinic unit cell if one exists; nothing otherwise\nwarn: Warn\nReturned lattice constants are not guaranteed to be standardized.\n\nExceptions\n\nThrows an ErrorException if the triclinic unit cell defined by lattice_constants is not equivalent to a base-centered monoclinic unit cell.\n\n\n\n\n\n","category":"function"},{"location":"functions/#XtallographyUtils.is_triclinic_type_I_cell","page":"Functions","title":"XtallographyUtils.is_triclinic_type_I_cell","text":"is_triclinic_type_I_cell(lattice_constants::TriclinicLatticeConstants) -> Bool\n\nDetermine whether the unit cell defined by lattice_constants is a Type I or Type II cell.\n\nA triclinic unit cell is Type I if the product of the dot products of all pairs of basis vectors for unit cell is positive:\n\n(veca cdot vecb)(vecb cdot vecc)(vecc cdot veca)  0\n\nOtherwise, the triclinic unit cell is Type II.\n\nReturn values\n\ntrue if lattice_constants defines a Type I cell; false if lattice_constants defines a Type II cell.\n\n\n\n\n\n","category":"function"},{"location":"functions/#XtallographyUtils.satisfies_triclinic_angle_constraints","page":"Functions","title":"XtallographyUtils.satisfies_triclinic_angle_constraints","text":"satisfies_triclinic_angle_constraints(α::Real, β::Real, γ::Real) -> Bool\n\nDetermine whether α, β, and γ satisfy the angle constraints for triclinic lattices:\n\n0   α + β + γ  2π\n0   α + β - γ  2π\n0   α - β + γ  2π\n0  -α + β + γ  2π\n\nReturn values\n\ntrue if (α, β, γ) form a valid triple of angles for a triclinic unit cell; false otherwise\n\nExamples\n\njulia> satisfies_triclinic_angle_constraints(π/4, π/5, π/6)\ntrue\njulia> satisfies_triclinic_angle_constraints(3π/4, 4π/5, 5π/6)\nfalse\n\n\n\n\n\n","category":"function"},{"location":"functions/#Monoclinic-Systems","page":"Functions","title":"Monoclinic Systems","text":"","category":"section"},{"location":"functions/","page":"Functions","title":"Functions","text":"convert_to_base_centering\nconvert_to_body_centering","category":"page"},{"location":"functions/#XtallographyUtils.convert_to_base_centering","page":"Functions","title":"XtallographyUtils.convert_to_base_centering","text":"convert_to_base_centering(\n    lattice_constants::MonoclinicLatticeConstants\n) -> MonoclinicLatticeConstants\n\nConvert a body-centered monoclinic unit cell to base-centered monoclinic unit cell.\n\nReturn values\n\nlattice constants for equivalent base-centered unit cell\n\nExamples\n\njulia> lattice_constants = MonoclinicLatticeConstants(1.0, 2.0, 3.0, 3π / 5);\n\njulia> base_centered_lattice_constants = convert_to_base_centering(lattice_constants);\n\njulia> base_centered_lattice_constants.a ≈ 2.8541019662496847\ntrue\n\njulia> base_centered_lattice_constants.b ≈ 2\ntrue\n\njulia> base_centered_lattice_constants.c ≈ 1\ntrue\n\njulia> base_centered_lattice_constants.β ≈ 1.5963584695539381\ntrue\n\n\n\n\n\n","category":"function"},{"location":"functions/#XtallographyUtils.convert_to_body_centering","page":"Functions","title":"XtallographyUtils.convert_to_body_centering","text":"convert_to_body_centering(\n    lattice_constants::MonoclinicLatticeConstants\n) -> MonoclinicLatticeConstants\n\nConvert a base-centered monoclinic unit cell to body-centered monoclinic unit cell.\n\nReturn values\n\nlattice constants for equivalent body-centered unit cell\n\nExamples\n\njulia> lattice_constants = MonoclinicLatticeConstants(1.0, 2.0, 3.0, 3π / 5);\n\njulia> body_centered_lattice_constants = convert_to_body_centering(lattice_constants);\n\njulia> body_centered_lattice_constants.a ≈ 2.8541019662496847\ntrue\n\njulia> body_centered_lattice_constants.b ≈ 2\ntrue\n\njulia> body_centered_lattice_constants.c ≈ 3\ntrue\n\njulia> body_centered_lattice_constants.β ≈ 2.8018712454717734\ntrue\n\n\n\n\n\n","category":"function"},{"location":"functions/","page":"Functions","title":"Functions","text":"","category":"page"},{"location":"functions/#Math-Functions","page":"Functions","title":"Math Functions","text":"","category":"section"},{"location":"functions/","page":"Functions","title":"Functions","text":"asin_\nacos_\nis_basis\nsurface_area(::Vector{<:Real}, ::Vector{<:Real}, ::Vector{<:Real})\nvolume(::Vector{<:Real}, ::Vector{<:Real}, ::Vector{<:Real})","category":"page"},{"location":"functions/#XtallographyUtils.asin_","page":"Functions","title":"XtallographyUtils.asin_","text":"asin_(x::Real; rtol::Real) -> Float64\n\nCompute the arcsin of x with a tolerance for values of x that are slightly outside of the mathematical domain [-1, 1].\n\nWhen the value of x is approximately equal to 1, asin_(x) returns π / 2; when the value of x is approximately equal to -1, asin_(x) returns -π / 2.\n\nReturn values\n\narcsin of x\n\nExamples\n\njulia> asin_(0.5) ≈ π / 6\ntrue\njulia> asin_(1 + eps(1.0)) == π / 2\ntrue\njulia> asin_(-1 - eps(1.0)) == -π / 2\ntrue\n\n\n\n\n\n","category":"function"},{"location":"functions/#XtallographyUtils.acos_","page":"Functions","title":"XtallographyUtils.acos_","text":"acos_(x::Real; rtol::Real) -> Float64\n\nCompute the arccos of x with a tolerance for values of x that are slightly outside of the mathematical domain [-1, 1].\n\nWhen the value of x is approximately equal to 1, acos_(x) returns 0; when the value of x is approximately equal to -1, acos_(x) returns π.\n\nReturn values\n\narccos of x\n\nExamples\n\njulia> acos_(0.5) ≈ π / 3\ntrue\njulia> acos_(1 + eps(1.0)) == 0\ntrue\njulia> acos_(-1 - eps(1.0)) == π\ntrue\n\n\n\n\n\n","category":"function"},{"location":"functions/#XtallographyUtils.is_basis","page":"Functions","title":"XtallographyUtils.is_basis","text":"is_basis(v1::Vector{Real}, v2::Vector{Real}, v3::Vector{Real}) -> Bool\n\nDetermine if the vectors v1, v2, and v3 are a basis for a three-dimensional lattice (i.e., v1, v2, and v3 are linearly independent).\n\nReturn values\n\ntrue if v1, v2, and v3 are a basis; false otherwise\n\nExamples\n\njulia> is_basis([1, 0, 0], [1, 1, 0], [1, 0, 1])\ntrue\n\njulia> is_basis([1, 0, 0], [1, 1, 0], [1, -1, 0])\nfalse\n\n\n\n\n\n","category":"function"},{"location":"functions/#XtallographyUtils.surface_area-Tuple{Vector{<:Real}, Vector{<:Real}, Vector{<:Real}}","page":"Functions","title":"XtallographyUtils.surface_area","text":"surface_area(v1::Vector{Real}, v2::Vector{Real}, v3::Vector{Real}) -> Float64\n\nCompute the surface area of the parallelipiped defined by the vectors v1, v2, and v3.\n\nReturn values\n\nsurface area of the parallelipiped defined by v1, v2, and v3\n\nExamples\n\njulia> surface_area([1, 0, 0], [1, 1, 0], [1, 0, 1])\n7.464101615137754\n\n\n\n\n\n","category":"method"},{"location":"functions/#XtallographyUtils.volume-Tuple{Vector{<:Real}, Vector{<:Real}, Vector{<:Real}}","page":"Functions","title":"XtallographyUtils.volume","text":"volume(v1::Vector{Real}, v2::Vector{Real}, v3::Vector{Real}) -> Float64\n\nCompute the volume of the parallelipiped defined by the vectors v1, v2, and v3.\n\nReturn values\n\nvolume of the parallelipiped defined by v1, v2, and v3\n\nExamples\n\njulia> volume([1, 0, 0], [1, 1, 0], [1, 0, 2])\n2.0\n\n\n\n\n\n","category":"method"},{"location":"functions/","page":"Functions","title":"Functions","text":"","category":"page"},{"location":"examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Create a unit cell object.\njulia> unit_cell = UnitCell(OrthorhombicLatticeConstants(2, 3, 4), body_centered)\nUnitCell(OrthorhombicLatticeConstants(2.0, 3.0, 4.0), BodyCentered())\nCompute basic unit cell attributes.\njulia> volume(unit_cell)\n24.0\n\njulia> surface_area(unit_cell)\n52.0\n\njulia> basis(unit_cell)\n([2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0])\nStandardize the lattice constants for a unit cell to be consistent with IUCr conventions.\njulia> unit_cell = UnitCell(OrthorhombicLatticeConstants(4, 2, 3), primitive)\nUnitCell(OrthorhombicLatticeConstants(4.0, 2.0, 3.0), Primitive())\n\njulia> standardize(unit_cell)\nUnitCell(OrthorhombicLatticeConstants(2.0, 3.0, 4.0), Primitive())\nCompute the Delaunay reduced cell for a unit cell.\njulia> unit_cell = UnitCell(CubicLatticeConstants(4), face_centered)\nUnitCell(CubicLatticeConstants(4.0), FaceCentered())\n\njulia> r_cell = reduced_cell(unit_cell)\nUnitCell(TriclinicLatticeConstants(2.8284271247461903, 2.8284271247461903, 2.8284271247461903, 1.0471975511965974, 1.0471975511965974, 1.5707963267948966), Primitive())\nCompute the conventional cell for a unit cell.\njulia> conventional_cell(r_cell)\nUnitCell(CubicLatticeConstants(3.9999999999999982), FaceCentered())","category":"page"},{"location":"docs-index/#Index","page":"Index","title":"Index","text":"","category":"section"},{"location":"docs-index/","page":"Index","title":"Index","text":"","category":"page"},{"location":"#XtallographyUtils","page":"Home","title":"XtallographyUtils","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"XtallographyUtils defines a collection of basic types and miscellaneous utility functions that support crystallography computations. The functionality provided by this package are not intended to be comprehensive. Rather, functionality is added on an as-needed basis to support research projects.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Currently, the Xtallography Utilities provides support for:","category":"page"},{"location":"","page":"Home","title":"Home","text":"types for defining Bravais lattice types,\nbasic unit cell computations (e.g., basis, volume, surface area),\nstandardization of lattice constants for unit cells, and\nconversions between equivalent unit cells for a lattice.","category":"page"},{"location":"#Getting-Started","page":"Home","title":"Getting Started","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Add the Velexi Julia package registry.\njulia>  # Press ']' to enter the Pkg REPL mode.\npkg> registry add https://github.com/velexi-research/JuliaRegistry.git\nnote: XtallographyUtils is registered with a local Julia package registry\nXtallographyUtils is registered with the Velexi Julia package registry (not the General Julia package registry), so the Pkg REPL will be able to find XtallographyUtils only if the Velexi Julia package registry has been added to your Julia installation. For more information about local registries for Julia packages, LocalRegistry.jl.\ntip: Only needed once\nThis step only needs to be performed once per Julia installation.\nInstall the XtallographyUtils package via the Pkg REPL. That's it!\njulia>  # Press ']' to enter the Pkg REPL mode.\npkg> add XtallographyUtils","category":"page"},{"location":"#Related-Packages","page":"Home","title":"Related Packages","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"There are a couple of crystallography packages in the Julia ecosystem that provide support for various crystall and lattice computations.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Crystallography.jl and related packages\nSpglib","category":"page"},{"location":"types/","page":"Types","title":"Types","text":"CurrentModule = XtallographyUtils","category":"page"},{"location":"types/#Types","page":"Types","title":"Types","text":"","category":"section"},{"location":"types/","page":"Types","title":"Types","text":"Lattice Types","category":"page"},{"location":"types/","page":"Types","title":"Types","text":"LatticeSystem\nCentering","category":"page"},{"location":"types/","page":"Types","title":"Types","text":"Unit Cell Types","category":"page"},{"location":"types/","page":"Types","title":"Types","text":"LatticeConstants\nUnitCell","category":"page"},{"location":"types/","page":"Types","title":"Types","text":"","category":"page"},{"location":"types/#Lattice-Types","page":"Types","title":"Lattice Types","text":"","category":"section"},{"location":"types/","page":"Types","title":"Types","text":"LatticeSystem\nCentering","category":"page"},{"location":"types/#XtallographyUtils.LatticeSystem","page":"Types","title":"XtallographyUtils.LatticeSystem","text":"LatticeSystem\n\nSupertype for the seven lattice systems in 3D\n\nSubtypes\n\nTriclinic, Monoclinic, Orthorhombic, Hexagonal, Rhombohedral, Tetragonal, Cubic\n\n\n\n\n\n","category":"type"},{"location":"types/#XtallographyUtils.Centering","page":"Types","title":"XtallographyUtils.Centering","text":"Centering\n\nSupertype for the four lattice centerings in 3D\n\nSubtypes\n\nPrimitive, BaseCentered, BodyCentered, FaceCentered\n\n\n\n\n\n","category":"type"},{"location":"types/#Concrete-Types","page":"Types","title":"Concrete Types","text":"","category":"section"},{"location":"types/#Lattice-Systems","page":"Types","title":"Lattice Systems","text":"","category":"section"},{"location":"types/","page":"Types","title":"Types","text":"All subtypes of LatticeSystem are singleton types. For convenience, a singleton instance is defined for each lattice system type.","category":"page"},{"location":"types/","page":"Types","title":"Types","text":"Triclinic\nMonoclinic\nOrthorhombic\nHexagonal\nRhombohedral\nTetragonal\nCubic","category":"page"},{"location":"types/#XtallographyUtils.Triclinic","page":"Types","title":"XtallographyUtils.Triclinic","text":"Triclinic\n\nType representing the triclinic lattice system that is the type of triclinic\n\nSupertype: LatticeSystem\n\n\n\n\n\n","category":"type"},{"location":"types/#XtallographyUtils.Monoclinic","page":"Types","title":"XtallographyUtils.Monoclinic","text":"Monoclinic\n\nType representing the monoclinic lattice system that is the type of monoclinic\n\nSupertype: LatticeSystem\n\n\n\n\n\n","category":"type"},{"location":"types/#XtallographyUtils.Orthorhombic","page":"Types","title":"XtallographyUtils.Orthorhombic","text":"Orthorhombic\n\nType representing the orthorhombic lattice system that is the type of orthorhombic\n\nSupertype: LatticeSystem\n\n\n\n\n\n","category":"type"},{"location":"types/#XtallographyUtils.Hexagonal","page":"Types","title":"XtallographyUtils.Hexagonal","text":"Hexagonal\n\nType representing the hexagonal lattice system that is the type of hexagonal\n\nSupertype: LatticeSystem\n\n\n\n\n\n","category":"type"},{"location":"types/#XtallographyUtils.Rhombohedral","page":"Types","title":"XtallographyUtils.Rhombohedral","text":"Rhombohedral\n\nType representing the rhombohedral lattice system that is the type of rhombohedral\n\nSupertype: LatticeSystem\n\n\n\n\n\n","category":"type"},{"location":"types/#XtallographyUtils.Tetragonal","page":"Types","title":"XtallographyUtils.Tetragonal","text":"Tetragonal\n\nType representing the tetragonal lattice system that is the type of tetragonal\n\nSupertype: LatticeSystem\n\n\n\n\n\n","category":"type"},{"location":"types/#XtallographyUtils.Cubic","page":"Types","title":"XtallographyUtils.Cubic","text":"Cubic\n\nType representing the cubic lattice system that is the type of cubic\n\nSupertype: LatticeSystem\n\n\n\n\n\n","category":"type"},{"location":"types/#Centerings","page":"Types","title":"Centerings","text":"","category":"section"},{"location":"types/","page":"Types","title":"Types","text":"All subtypes of Centering are singleton types. For convenience, a singleton instance is defined for each centering type.","category":"page"},{"location":"types/","page":"Types","title":"Types","text":"Primitive\nBaseCentered\nBodyCentered\nFaceCentered","category":"page"},{"location":"types/#XtallographyUtils.Primitive","page":"Types","title":"XtallographyUtils.Primitive","text":"Primitive\n\nType representing no centering that is the type of primitive\n\nSupertype: Centering\n\n\n\n\n\n","category":"type"},{"location":"types/#XtallographyUtils.BaseCentered","page":"Types","title":"XtallographyUtils.BaseCentered","text":"BaseCentered\n\nType representing base centering that is the type of base_centered\n\nnote: Note\nBy convention, base-centering is on the C-face of the unit cell.\n\nSupertype: Centering\n\n\n\n\n\n","category":"type"},{"location":"types/#XtallographyUtils.BodyCentered","page":"Types","title":"XtallographyUtils.BodyCentered","text":"BodyCentered\n\nType representing body centering that is the type of body_centered\n\nSupertype: Centering\n\n\n\n\n\n","category":"type"},{"location":"types/#XtallographyUtils.FaceCentered","page":"Types","title":"XtallographyUtils.FaceCentered","text":"FaceCentered\n\nType representing face centering that is the type of face_centered\n\nSupertype: Centering\n\n\n\n\n\n","category":"type"},{"location":"types/#Constants","page":"Types","title":"Constants","text":"","category":"section"},{"location":"types/#Lattice-Systems-2","page":"Types","title":"Lattice Systems","text":"","category":"section"},{"location":"types/","page":"Types","title":"Types","text":"triclinic\nmonoclinic\northorhombic\nhexagonal\nrhombohedral\ntetragonal\ncubic","category":"page"},{"location":"types/#XtallographyUtils.triclinic","page":"Types","title":"XtallographyUtils.triclinic","text":"triclinic\n\nThe singleton instance of type Triclinic\n\n\n\n\n\n","category":"constant"},{"location":"types/#XtallographyUtils.monoclinic","page":"Types","title":"XtallographyUtils.monoclinic","text":"monoclinic\n\nThe singleton instance of type Monoclinic\n\n\n\n\n\n","category":"constant"},{"location":"types/#XtallographyUtils.orthorhombic","page":"Types","title":"XtallographyUtils.orthorhombic","text":"orthorhombic\n\nThe singleton instance of type Orthorhombic\n\n\n\n\n\n","category":"constant"},{"location":"types/#XtallographyUtils.hexagonal","page":"Types","title":"XtallographyUtils.hexagonal","text":"hexagonal\n\nThe singleton instance of type Hexagonal\n\n\n\n\n\n","category":"constant"},{"location":"types/#XtallographyUtils.rhombohedral","page":"Types","title":"XtallographyUtils.rhombohedral","text":"rhombohedral\n\nThe singleton instance of type Rhombohedral\n\n\n\n\n\n","category":"constant"},{"location":"types/#XtallographyUtils.tetragonal","page":"Types","title":"XtallographyUtils.tetragonal","text":"tetragonal\n\nThe singleton instance of type Tetragonal\n\n\n\n\n\n","category":"constant"},{"location":"types/#XtallographyUtils.cubic","page":"Types","title":"XtallographyUtils.cubic","text":"cubic\n\nThe singleton instance of type Cubic\n\n\n\n\n\n","category":"constant"},{"location":"types/#Centerings-2","page":"Types","title":"Centerings","text":"","category":"section"},{"location":"types/","page":"Types","title":"Types","text":"primitive\nbase_centered\nbody_centered\nface_centered","category":"page"},{"location":"types/#XtallographyUtils.primitive","page":"Types","title":"XtallographyUtils.primitive","text":"primitive\n\nThe singleton instance of type Primitive\n\n\n\n\n\n","category":"constant"},{"location":"types/#XtallographyUtils.base_centered","page":"Types","title":"XtallographyUtils.base_centered","text":"base_centered\n\nThe singleton instance of type BaseCentered\n\n\n\n\n\n","category":"constant"},{"location":"types/#XtallographyUtils.body_centered","page":"Types","title":"XtallographyUtils.body_centered","text":"body_centered\n\nThe singleton instance of type BodyCentered\n\n\n\n\n\n","category":"constant"},{"location":"types/#XtallographyUtils.face_centered","page":"Types","title":"XtallographyUtils.face_centered","text":"face_centered\n\nThe singleton instance of type FaceCentered\n\n\n\n\n\n","category":"constant"},{"location":"types/#Other-Constants","page":"Types","title":"Other Constants","text":"","category":"section"},{"location":"types/","page":"Types","title":"Types","text":"BRAVAIS_LATTICES","category":"page"},{"location":"types/#XtallographyUtils.BRAVAIS_LATTICES","page":"Types","title":"XtallographyUtils.BRAVAIS_LATTICES","text":"BRAVAIS_LATTICES\n\nList of valid Bravais lattices.\n\n\n\n\n\n","category":"constant"},{"location":"types/","page":"Types","title":"Types","text":"","category":"page"},{"location":"types/#Unit-Cell-Types","page":"Types","title":"Unit Cell Types","text":"","category":"section"},{"location":"types/","page":"Types","title":"Types","text":"LatticeConstants\nUnitCell","category":"page"},{"location":"types/#XtallographyUtils.LatticeConstants","page":"Types","title":"XtallographyUtils.LatticeConstants","text":"LatticeConstants\n\nSupertype for lattice constants for the seven lattice systems in 3D\n\nSubtypes\n\nTriclinicLatticeConstants, MonoclinicLatticeConstants, OrthorhombicLatticeConstants, TetragonalLatticeConstants, RhombohedralLatticeConstants, HexagonalLatticeConstants, CubicLatticeConstants\n\n\n\n\n\n","category":"type"},{"location":"types/#XtallographyUtils.UnitCell","page":"Types","title":"XtallographyUtils.UnitCell","text":"UnitCell\n\nUnit cell for a lattice\n\nFields\n\nlattice_constants: lattice constants of unit cell\ncentering: centering of unit cell\n\n\n\n\n\n","category":"type"},{"location":"types/#Lattice-Constants","page":"Types","title":"Lattice Constants","text":"","category":"section"},{"location":"types/","page":"Types","title":"Types","text":"TriclinicLatticeConstants\nMonoclinicLatticeConstants\nOrthorhombicLatticeConstants\nHexagonalLatticeConstants\nRhombohedralLatticeConstants\nTetragonalLatticeConstants\nCubicLatticeConstants","category":"page"},{"location":"types/#XtallographyUtils.TriclinicLatticeConstants","page":"Types","title":"XtallographyUtils.TriclinicLatticeConstants","text":"TriclinicLatticeConstants\n\nLattice constants for a triclinic unit cell\n\nFields\n\na, b, c: lengths of the edges of the unit cell\nα, β, γ: angles between edges of the unit cell in the planes of the faces of the unit cell\n\nnote: Note\nThe constraints that valid triclinic unit cell angles must satisfy are not enforced by the constructor. It is acceptable to construct TriclinicLatticeConstants with invalid values for α, β, and γ.\n\nSupertype: LatticeConstants\n\n\n\n\n\n","category":"type"},{"location":"types/#XtallographyUtils.MonoclinicLatticeConstants","page":"Types","title":"XtallographyUtils.MonoclinicLatticeConstants","text":"MonoclinicLatticeConstants\n\nLattice constants for a monoclinic unit cell\n\nFields\n\na, b, c: lengths of the edges of the unit cell\nβ: angle between edges of the unit cell in the plane of the face of the unit cell where the edges are not orthogonal\n\nSupertype: LatticeConstants\n\n\n\n\n\n","category":"type"},{"location":"types/#XtallographyUtils.OrthorhombicLatticeConstants","page":"Types","title":"XtallographyUtils.OrthorhombicLatticeConstants","text":"OrthorhombicLatticeConstants\n\nLattice constants for an orthorhombic unit cell\n\nFields\n\na, b, c: lengths of the edges of the unit cell\n\nnote: Note\nBy convention, edge lengths of the unit cell are ordered so that a <= b <= c.\n\nSupertype: LatticeConstants\n\n\n\n\n\n","category":"type"},{"location":"types/#XtallographyUtils.HexagonalLatticeConstants","page":"Types","title":"XtallographyUtils.HexagonalLatticeConstants","text":"HexagonalLatticeConstants\n\nLattice constants for a hexagonal unit cell\n\nFields\n\na, c: lengths of the edges of the unit cell\n\nSupertype: LatticeConstants\n\n\n\n\n\n","category":"type"},{"location":"types/#XtallographyUtils.RhombohedralLatticeConstants","page":"Types","title":"XtallographyUtils.RhombohedralLatticeConstants","text":"RhombohedralLatticeConstants\n\nLattice constants for a rhombohedral unit cell\n\nFields\n\na: length of the edge of the unit cell\nα: angle between edges of the unit cell in the plane of the faces of the unit cell\n\nSupertype: LatticeConstants\n\n\n\n\n\n","category":"type"},{"location":"types/#XtallographyUtils.TetragonalLatticeConstants","page":"Types","title":"XtallographyUtils.TetragonalLatticeConstants","text":"TetragonalLatticeConstants\n\nLattice constants for a tetragonal unit cell\n\nFields\n\na, c: lengths of the edges of the unit cell\n\nSupertype: LatticeConstants\n\n\n\n\n\n","category":"type"},{"location":"types/#XtallographyUtils.CubicLatticeConstants","page":"Types","title":"XtallographyUtils.CubicLatticeConstants","text":"CubicLatticeConstants\n\nLattice constants for a cubic unit cell\n\nFields\n\na: length of the edge of the unit cell\n\nSupertype: LatticeConstants\n\n\n\n\n\n","category":"type"},{"location":"types/","page":"Types","title":"Types","text":"","category":"page"}]
}
