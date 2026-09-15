nextflow.enable.dsl=2

process PREPARE_RECEPTORS {
        tag "Fixing PDB: ${pdb_id}"
        publishDir "results/receptors_fixed", mode: 'copy'
        container 'quay.io/biocontainers/pdbfixer:1.9--pyh5e36f58_0'

        input:
        val pdb_id

        output:
        tuple val(pdb_id), path("${pdb_id}_fixed.pdb"), emit: fixed_bundle

        script:
        """      
        pdbfixer --pdbid=${pdb_id} --output=${pdb_id}_fixed.pdb --add-atoms=heavy --replace-nonstandard
        """
}

process GENERATE_TOPOLOGY {
        tag "Building Topology: ${pdb_id}"
        publishDir "results/topology_generated", mode: 'copy'
        container 'ghcr.io/mscbioinformatics/vmd:1.9.4'

        input:
        tuple val(pdb_id), path(fixed_pdb)
        path toppar_dir

        output:
        tuple val(pdb_id), path("${pdb_id}_solv_ion.pdb"), path("${pdb_id}_solv_ion.psf"), emit: ready_for_namd

        script:
        """
        cat << 'EOF' > generate_psf.tcl
        package require psfgen
        topology ${toppar_dir}/top_all36_prot.rtf
        
        segment PROT {
                 pdb ${fixed_pdb}
        }
        coordpdb ${fixed_pdb} PROT
        guesscoord
        writepdb "APO_prot.pdb"
        writepsf "APO_prot.psf"
        package require solvate
        solvate APO_prot.psf APO_prot.pdb -t 10 -o solvated
    
# 1. Solvatación con colchón (12 Ä)

        package require autoionize
        solvate APO_prot.psf APO_prot.pdb -t 12 -o solvated

# 2. Ionización y Neutralización elecostática (0.15 M NaCl)

        autoionize -psf solvated.psf -pdb solvated.pdb -sc 0.15 -neutral -o ${pdb_id}_solv_ion
        exit
        EOF

        vmd -dispdev text -e generate_psf.tcl
        """
}

process WRITE_NAMD_CONFIGS {
        tag "NAMD Configs/PME Generation: ${pdb_id}"
        publishDir "results/namd_inputs", mode: 'copy'
        container 'python:3.10-slim'
    
        input:
        tuple val(pdb_id), path(solv_pdb), path(solv_psf), val(replica)

        output:
        tuple val(pdb_id), path(solv_pdb), path(solv_psf), path("${pdb_id}_r${replica}_min_nvt.conf"), path("${pdb_id}_r${replica}_npt.conf"), path("${pdb_id}_r${replica}_equi.conf"), path("${pdb_id}_r${replica}_prod.conf"), emit: Configs

        script:

        """
        cat << 'EOF' > calculate_pme.py
        import numpy as np

        x, y, z = [], [], []
        with open("${solv_pdb}", "r") as f:
                for line in f:
                        if line.startswith("ATOM") or line.startswith("HETATM"):
                                x.append(float(line[30:38]))
                                y.append(float(line[38:46]))
                                z.append(float(line[46:54]))

        caja_x, caja_y, caja_z = max(x)-min(x), max(y)-min(y), max(z)-min(z)

        def proximo_primo_pme(n):
                n = int(np.ceil(n))
                while True:
                        m = n
                        for p in:
                                while m % p == 0: m //= p
                if m == 1: return n
                n += 1

        pme_x, pme_y, pme_z = proximo_primo_pme(caja_x), proximo_primo_pme(caja_y), proximo_primo_pme(caja_z)

# 1. Parámetros físicos genelares en común para NAMD

        def escribir_base(f):
                f.write("structure          ${solv_psf}\\n")
                f.write("coordinates        ${solv_pdb}\\n")
                f.write("cutoff             12.0\\n")
                f.write("switching          on\\n")
                f.write("switchdist         10.0\\n")
                f.write("pairlistdist       14.0\\n")
                f.write("PME                yes\\n")
                f.write(f"PMEGridSizeX       {pme_x}\\n")
                f.write(f"PMEGridSizeY       {pme_y}\\n")
                f.write(f"PMEGridSizeZ       {pme_z}\\n")
                f.write("timestep           2.0\\n")
                f.write("rigidBonds         all\\n")
                f.write("nonbondedFreq      1\\n")
                f.write("fullElectFrequency 2\\n")
                f.write("stepspercycle      10\\n")
                f.write(f"seed              {12345 + int(replica)}\\n")

 # 2. Fase de MINIMIZACION

        with open("${pdb_id}_r${replica}_min_nvt.conf", "w") as f:
                escribir_base(f)
                f.write("temperature        310\\n")
                f.write(f"cellBasisVector1   {caja_x:.2f} 0.0 0.0\\n")
                f.write(f"cellBasisVector2   0.0 {caja_y:.2f} 0.0\\n")
                f.write(f"cellBasisVector3   0.0 0.0 {caja_z:.2f}\\n")
                f.write(f"cellOrigin         {np.mean(x):.2f} {np.mean(y):.2f} {np.mean(z):.2f}\\n")
                f.write("langevin           on\\nlangevinDamping    1.0\\nlangevinTemp       310\\n")
                f.write("minimize           1000\\n")
                f.write("reinitvels         310\\n")
                f.write("run                50000\\n") # 100 ps NVT

 # 3. Fase de ANNEALING (NVT)
    
        with open("${pdb_id}_r${replica}_equi.conf", "w") as f:
                escribir_base(f)
                f.write("binCoordinates     ${pdb_id}_r${replica}_min_equi.coor\\n")
                f.write("binVelocities      ${pdb_id}_r${replica}_min_equi.vel\\n")
                f.write("extendedSystem     ${pdb_id}_r${replica}_min_equi.xsc\\n")
                f.write("langevin           on\\nlangevinDamping    1.0\\nlangevinTemp       310\\n")
                f.write("LangevinPiston     on\\nLangevinPistonTarget 1.01325\\nLangevinPistonPeriod 200\\nLangevinPistonDecay  100\\nLangevinPi>
                f.write("useGroupPressure   yes\\nuseFlexibleCell   no\\n") # no para sistemas sin membrana
                f.write("run                100000\\n") # 200 ps NPT

 # 4. Fase de EQUILIBRATION
        
        with open("${pdb_id}_npt.conf", "w") as f:
                escribir_base(f)
                f.write("binCoordinates     ${pdb_id}_min_nvt.coor\\n")
                f.write("binVelocities      ${pdb_id}_min_nvt.vel\\n")
                f.write("extendedSystem     ${pdb_id}_min_nvt.xsc\\n")
                f.write("langevin           on\\nlangevinDamping    1.0\\nlangevinTemp       310\\n")
                f.write("LangevinPiston     on\\nLangevinPistonTarget 1.01325\\nLangevinPistonPeriod 200\\nLangevinPistonDecay  100\\nLangevinPi>
                f.write("useGroupPressure   yes\\nuseFlexibleCell   no\\n") # no para sistemas sin membrana
                f.write("run                100000\\n") # 200 ps NPT

# 5. Fase de PRODUCCION (100 ns)
       
        with open("${pdb_id}_prod.conf", "w") as f:
                escribir_base(f)
        
        # Continuidad desde NPT con checkpoints binarios
        
                
                f.write(f"seed               {12345 + replica}\\n")
                f.write("binCoordinates     ${pdb_id}_npt.restart.coor\\n")
                f.write("binVelocities      ${pdb_id}_npt.restart.vel\\n")
                f.write("extendedSystem     ${pdb_id}_npt.restart.xsc\\n")
                f.write("firstTimestep      0\\n")

        # Termostato (Langevin)

                f.write("langevin           on\\n") 
                f.write("langevinDamping    1.0\\n")
                f.write("langevinTemp       310\\n")
                f.write("langevinHydrogen   off\\n")
        
        # Barostato (Langevin Piston, NPT)

                f.write("langevinPiston          on\\n") 
                f.write("langevinPistonTarget    1.0\\n")
                f.write("langevinPistonPeriod    200.0\\n")
                f.write("langevinPistonDecay     100.0\\n")
                f.write("langevinPistonTemp      310\\n")
                f.write("useGroupPressure        yes\\n")
                f.write("useFlexibleCell         no\\n")
                f.write("useConstantArea         no\\n")
               
        # Aceleración por GPU

                f.write("gpuResident"    on\\n")
                f.write("CUDAMode"       on\\n")

        # Integrador

                f.write("timestep           2.0\\n")
                f.write("fullElectFrequency 2\\n")
                f.write("nonbondedFreq      1\\n")
                f.write("stepsPerCycle      200\\n")
                f.write("pairlistsPerCycle  2\\n")
               
        # Reproducibilidad

        f.write("seed                       12345\\n") # cambia por réplica

        # Output: trayectoria

                f.write("outputName         ${pdb_id}_prod\\n")
                f.write("dcdfile            ${pdb_id}_prod.dcd\\n")
                f.write("dcdFreq            10\\n")
                f.write("XSTFreq            10\\n")
                
                f.write("binaryOutput       no\\n")
    
        # Output: checkpoints

                f.write("restartFreq        10\\n")
                f.write("restartName        ${pdb_id}_prod.restart\\n")
                f.write("binaryRestart      yes\\n")

        # Output: log para evitar la limitacion de I/O

                f.write("outputEnergies     40\\n")
                f.write("outputTiming       40\\n")
                f.write("outputPressure     40\\n")
                
        # Duración de la dinámica: 100 ns
               
                f.write("run                50000000\\n")

        EOF
        python3 calculate_pme.py
        """
}

process RUN_NAMD_PIPELINE {
        tag "Running ${etapa_actual} de NAMD: ${pdb_id}"
        publishDir "results/namd_outputs", mode: 'copy'
        container 'nvcr.io/hpc/namd:3.0'

        input:
        tuple val(pdb_id), path(solv_pdb), path(solv_psf), path(min_conf), path(npt_conf), path(equi_conf), path(prod_conf)

        output:
        path "${pdb_id}_prod.dcd", emit: trajectory
        path "${pdb_id}_prod.log", emit: log

        script:

         """
 # 1. Ejecucion secuencial en GPU nativa

        namd3 +p2 +setcpuaffinity +devices 0 +CUDASOA 1 ${min_conf} > ${pdb_id}_min_nvt.log
        namd3 +p2 +setcpuaffinity +devices 0 +CUDASOA 1 ${npt_conf} > ${pdb_id}_npt.log
        namd3 +p2 +setcpuaffinity +devices 0 +CUDASOA 1 ${equi_conf} > ${pdb_id}_equi.log
        namd3 +p2 +setcpuaffinity +devices 0 +CUDASOA 1 ${prod_conf} > ${pdb_id}_prod.log
        """
}

process RUN_BIO3D_PCA_KMEANS {
        tag "Running Bio3D for trayectory processing of: ${pdb_id}"
        publishDir "results/Bio3D_outputs", mode: 'copy'

        input:
        tuple val(pdb_id), path(solv_pdb), path(solv_psf), path(min_conf), path(npt_conf), path(equi_conf), path(prod_conf)

        output:
        path "${pdb_id}_medoid_cluster.*.pdb.", emit:

        script:

         """
        # 1. PCA, Paisaje de Energía libre y Clustering K-Means (medoides)

        library(bio3d) 
        
        # 1 Cargar Estructura de Referencia (PSF/PDB) y Trayectoria (DCD)
        
        psf_file <- "data/fimh_apo.psf"
        dcd_file <- "data/fimh_apo_prod.dcd"
        
        pdb <- read.pdb(psf_file)
        dcd <- read.dcd(dcd_file)
        
        # 2. Seleccionar Carbonos Alfa (CA) para eliminar ruido térmico
        
        ca.inds <- atom.select(pdb, "calpha")
        
        # 3. Superposición/Alineación Estructural sobre el primer frame
        
        xyz <- fit.xyz(fixed = pdb$xyz, mobile = dcd,
                fixed.inds = ca.inds$xyz, mobile.inds = ca.inds$xyz)
        
        # 4. Análisis de Componentes Principales (PCA)
        
        pc <- pca.xyz(xyz[, ca.inds$xyz])
        
        # Gráfico de Varianza Acumulada (Scree Plot)
        
        pdf("results/pca_scree_plot.pdf")
        plot(pc, pc.axes=1:2, main="PCA - ${pdb_id} APO Trajectory")
        dev.off()
        
        # 5. Cálculo del Paisaje de Energía Libre (FEL) sobre PC1 y PC2
                # G = -kB * T * ln(P(PC1, PC2))
        pc1 <- pc$z[, 1]
        pc2 <- pc$z[, 2]
        
        # 6. Clustering K-means sobre el Espacio de Componentes Principales
        
        set.seed(42)
        k_clusters <- 3
        km <- kmeans(cbind(pc1, pc2), centers = k_clusters)
        
        # 7. Identificación del MEDOIDE (Frame real más cercano al centro de cada clúster)
        
        medoid_frames <- c()
        for(i in 1:k_clusters) {
        cluster_indices <- which(km$cluster == i)
        cluster_coords  <- cbind(pc1, pc2)[cluster_indices, ]
        center <- km$centers[i, ]
        
        # Distancia euclidiana al centroide del clúster
        
        distances <- sqrt((cluster_coords[,1] - center[1])^2 + (cluster_coords[,2] - center[2])^2) medoid_frames[i] <- cluster_indices[which.min(distances)]
        }
        
        cat("Frames Medoides Seleccionados para Ensemble Docking:", medoid_frames, "\\n")
        
        # 8. Exportar los PDBs de los Medoides Reales para Virtual Screening
        
        for(i in 1:length(medoid_frames)) {
        frame_num <- medoid_frames[i]
        write.pdb(pdb = pdb, xyz = xyz[frame_num, ],
                file = paste0("results/medoid_cluster_", i, "_frame_", frame_num, ".pdb"))
        }
        
        cat("¡Medoides exportados exitosamente en 'results/Medoides'!\\n")

        """
}

workflow {
        main:
        raw_pdbs = Channel.of('9AT9', '9HW9')
        replicas_HOLO = Channel.of(1, 2, 3)

        PREPARE_RECEPTORS(raw_pdbs)
        GENERATE_TOPOLOGY(PREPARE_RECEPTORS.out.fixed_bundle)
        topology_with_replicas = GENERATE_TOPOLOGY.out.ready_for_namd.combine(replicas_HOLO)
        WRITE_NAMD_CONFIG(topology_with_replicas)
        RUN_NAMD_SIMULATION(WRITE_NAMD_CONFIG.out.configs)
}
