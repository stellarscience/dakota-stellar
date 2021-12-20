CreateTestGroup(
    MODELING_CORE

    WorkPieceTests.cpp
    #WorkPieceDerivativesTests.cpp
    WorkGraphTests.cpp
    WorkGraphPieceTests.cpp
    ModGraphPieceTests.cpp
    ModPieceTests.cpp
    ConstantPieceTests.cpp
    DerivativePieceTests.cpp
    IdentityPieceTests.cpp
    SplitVectorTests.cpp
    ScaleVectorTests.cpp
    CombineVectorsTests.cpp
    SumPieceTests.cpp
    MultiLogisticTest.cpp
)

CreateTestGroup(
  MODELING_STAN
  CwiseOperatorTests/UnaryOperatorTest.cpp
)

CreateTestGroup(
    MODELING_FLANN

    FlannTest/FlannCacheTests.cpp
)

CreateTestGroup(
    MODELING_LINEARALGEBRA

    LinearAlgebra/LinearOperatorTests.cpp
    LinearAlgebra/LyapunovSolverTests.cpp
    LinearAlgebra/AffineOperatorTests.cpp
    LinearAlgebra/LOBPCGTests.cpp
    LinearAlgebra/StochasticEigenSolverTests.cpp
    LinearAlgebra/SliceOperatorTests.cpp

    #LinearAlgebra/ScalarAlgebraTests.cpp
    #LinearAlgebra/EigenVectorAlgebraTests.cpp
    #LinearAlgebra/EigenMatrixAlgebraTests.cpp
    #LinearAlgebra/SundialsAlgebraTests.cpp
)

CreateTestGroup(
    MODELING_SUNDIALS_MODELS

    ODETests.cpp
    #RootfindingIVPTests.cpp
)

CreateTestGroup(
    MODELING_DISTRIBUTIONS

    DistributionsTests/DistributionTests.cpp
    DistributionsTests/InverseGammaTests.cpp
    DistributionsTests/UniformDistributionTests.cpp
    DistributionsTests/GaussianDistributionTests.cpp
    DistributionsTests/DensityProductTests.cpp
)
