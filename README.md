# Tooth-molar-ssm
Tooth-molar-ssm/main/scala/api/sampling/MixedProposalDistributions.scala
Tooth-molar-ssm/main/scala/api/sampling/evaluators/IndependentPointDistanceEvaluator.scala
Tooth-molar-ssm/main/scala/apps/teeth/registeration/

Tooth-molar-ssm/main/scala/apps/teeth/definition/changepointKernelForPartialMesh.scala
Tooth-molar-ssm/main/scala/apps/teeth/definition/changepointKernel.scala

 //load full samples data
    val (model, modelLms, sampleMeshes, sampleLMs) = LoadTestData.modelAndTarget("Full")
    //align meshes
    val Seq(alignedMesh,alignedLm) = alignShapes(modelLms,sampleMeshes,sampleLMs).align()

    //registeration

    //--non rigid icp
    //--grad opt
    //--mcmc
    //build pca model
    //Augment pca model
    //load partial data
    //partial registeration
    //--non rigid icp
    //-- icp proposal mcmc
    //pca model
