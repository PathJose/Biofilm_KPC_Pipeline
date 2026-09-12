nextflow.enable.dsl=2

process PREPARE_RECEPTORS {
        tag "Fixing PDB: ${pdb_id}"
        publishDir "results/receptors_fixed", mode: 'copy'
        conda 'bioconda::pdbfixer=1.9'

        input:
        val pdb_id

        output:
        path "${pdb_id}_fixed.pdb", emit: fixed_pdb

        script:
        """
        pdbfixer --pdbid=${pdb_id} --output=${pdb_id}_fixed.pdb --add-atoms=heav
        """
}

workflow {
        main:
        raw_pdbs = Channel.of('9AT9', '9HW9')

        PREPARE_RECEPTOR(raw_pdbs)
}
