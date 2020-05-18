# AlpineBugs_PNAS2020
Computer Code used to analyze the data and produce results for: 

Muhlfeld, CC, Cline, TJ, Giersch, JJ, Florentine, C, Pietzsch, E., Jacobsen, D., Hotaling, S.
2020. 
Rare meltwater biodiversity persists despite widespread deglaciation. 
Proceedings of the National Academy of Sciences.

AlpineBugs_Regression_PNAS.R - Calculates regression analyses of Richness ~ Distance to source, elevation, glacier cover presented in Fig 2.

AlpineBugs_LDA_and_NoGlacierAnalysis_PNAS.R - Partitions communities based on abudance using Latent Dirichlet Allocation and calculates regressions for sites without glaciers. These analyses support figures and tables related to Fig. 3.

AlpineBugs_FutureProjections_PNAS.R - Uses weighted model predictions to project habitat availability for cold-water community under different temperature scenarios. Results presented in Fig. 4.

AlpineBugs_BetaPartitioning_PNAS.R - Partitions beta-diversity into nestedness and turnover components following Balsega 2010. Shows that communities differ among sites rather than nested component communities.

AlpineBugs_IndividaulRarefactionCurves_PNAS.R - Calculates rarefied 'sample' to determine whether missing rare species may influence regression of Rich ~ Glacier Cover, distance, elev. Rarefication did not change relationships presented in Fig. 2.
