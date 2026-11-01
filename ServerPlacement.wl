(* ::Package:: *)

(*Creates a function that allows variables to be changed and accordingly changes the manipulate*)

Clear[ServerPlacementSimulation];

Options[ServerPlacementSimulation] = {"ServerRadius" -> 0.015, "HeatGeneration" -> 1000, "MaxTime" -> 600};

ServerPlacementSimulation[n_?IntegerQ, roomWidth_?NumericQ, roomHeight_
    ?NumericQ, plotTime_?NumericQ, opts : OptionsPattern[]] :=
    Module[{spec, serverRadius, minX, maxX, minY, maxY},
        serverRadius = OptionValue[ServerPlacementSimulation, "ServerRadius"];
        maxX = roomWidth;
        maxY = roomHeight;
        {minX, minY} = {0, 0};
        spec = Transpose[{RandomReal[{minX + serverRadius, maxX - serverRadius}, n], RandomReal[{minY + serverRadius, maxY - serverRadius}, n]}];
        ServerPlacementSimulation[spec, roomWidth, roomHeight, plotTime, opts]
    ];

(*ServerPlacementSimulation[spec:{_?NumericQ,_?NumericQ},roomWidth_?NumericQ,roomHeight_?NumericQ,plotTime_?NumericQ,opts:OptionsPattern[]]:=
ServerPlacementSimulation[{spec},roomWidth,roomHeight,plotTime,opts];*)

ServerPlacementSimulation[spec : {{_?NumericQ, _?NumericQ}..}, roomWidth_
    ?NumericQ, roomHeight_?NumericQ, plotTime_?NumericQ, opts : OptionsPattern[
    ]] :=
    Module[{minX = 0, minY = 0, maxX, maxY, serverRadius, heatGeneration,
         tend, \[CapitalOmega]2D, serverRegions, sourceRegion, rmf, vars, pars, pde, Tfun, 
        kAir, cpAir, \[Rho]Air},
        serverRadius = OptionValue[ServerPlacementSimulation, "ServerRadius"];
        heatGeneration = OptionValue[ServerPlacementSimulation, "HeatGeneration"];
        tend = OptionValue[ServerPlacementSimulation, "MaxTime"];
        maxX = roomWidth;
        maxY = roomHeight;
        kAir = QuantityMagnitude[ThermodynamicData["Air", "ThermalConductivity"
            ]];
        cpAir = QuantityMagnitude[ThermodynamicData["Air", "IsobaricHeatCapacity"
            ]];
        \[Rho]Air = QuantityMagnitude[ThermodynamicData["Air", "Density"]]
            ;
        vars = {T[t, x, y], t, {x, y}};
        \[CapitalOmega]2D = Rectangle[{minX, minY}, {maxX, maxY}];
        serverRegions = Rectangle[# - {serverRadius, serverRadius}, #
             + {serverRadius, serverRadius}]& /@ spec;
        sourceRegion = RegionUnion @@ serverRegions;
        rmf = RegionMember[sourceRegion];
        pars =
            <|
                "ThermalConductivity" -> kAir
                ,
                "MassDensity" -> \[Rho]Air
                ,
                "SpecificHeatCapacity" -> cpAir
                ,
                "HeatSource" ->
                    With[{Qg = heatGeneration, rmfLocal = rmf},
                        If[rmfLocal[{x, y}],
                            Qg
                            ,
                            0
                        ]
                    ]
            |>;
        pde = {HeatTransferPDEComponent[vars, pars] == NeumannValue[0,
             x == minX || x == maxX || y == minY || y == maxY], T[0, x, y] == 22}
            ;
        Tfun = NDSolveValue[pde, T, {t, 0, tend}, {x, y} \[Element] \[CapitalOmega]2D];
        <|"Spec" -> spec, "Solution" -> Tfun, "RoomHeight" -> maxY, "RoomWidth"
             -> maxX|>
    ];
    
Clear[SymmetricConfigurations];
SymmetricConfigurations[conf_:{{_,_}..},nRows_Integer,nColumns_Integer]:=
	SymmetricConfigurations[conf,nRows,nColumns]=
	 Module[{gridPoints,reflectedPoints,horizontallySymmetricalPoints,verticallySymmetricalPoints},
		gridPoints = Flatten[Table[{i,j},{i,nRows},{j,nColumns}],1];
		reflectedPoints={nRows-#[[1]]+1,#[[2]]}&/@conf;
		horizontallySymmetricalPoints = Select[gridPoints,MemberQ[reflectedPoints,#]&];
		reflectedPoints={#[[1]],nColumns-#[[2]]+1}&/@conf;
		verticallySymmetricalPoints = Select[gridPoints,MemberQ[reflectedPoints,#]&];
		{horizontallySymmetricalPoints,verticallySymmetricalPoints}
	];

Clear[SameCombinationQ];
SameCombinationQ[comb1:{{_?NumericQ,_?NumericQ}..},comb2:{{_?NumericQ,_?NumericQ}..},nRows_Integer,nColumns_Integer]:=
	comb1==comb2||MemberQ[Sort/@SymmetricConfigurations[comb1,nRows,nColumns],Sort@comb2]
	
Clear[UniqueOptimizationCandidates];
UniqueOptimizationCandidates[gridPoints_,k_,nRows_,nColumns_]:=
	Module[{validCombinations,comb2,comb3,comb4},validCombinations=Tuples[gridPoints,k];
		comb2 = Select[validCombinations,Length[Union[#]]==Length[#]&];
		comb3=Union[Sort/@comb2];
		comb4=Union[comb3,SameTest->(SameCombinationQ[#1,#2,nRows,nColumns]&)];
		comb4
	]
