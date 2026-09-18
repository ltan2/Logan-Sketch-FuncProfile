COPY (
  SELECT acc
  FROM sra_metadata
  WHERE UPPER(TRIM(assay_type))    = 'WGS'
    AND UPPER(TRIM(librarysource)) = 'METAGENOMIC'
    AND (
          -- tier 1: organism taxonomy says metagenome
          organism ILIKE '%metagenome%'

          -- tier 2: BioSample checklist is metagenome-specific
       OR array_to_string(biosamplemodel_sam, ',') ILIKE '%MIMAG%'
       OR array_to_string(biosamplemodel_sam, ',') ILIKE '%MIMS.me%'
       OR array_to_string(biosamplemodel_sam, ',') ILIKE '%MIUVIG%'
       OR array_to_string(biosamplemodel_sam, ',') ILIKE '%Metagenome or environmental%'

          -- tier 3: MIxS environmental fields
       OR (
            (    jattr ILIKE '%env_broad_scale%'
              OR jattr ILIKE '%env_local_scale%'
              OR jattr ILIKE '%env_medium%' )
          )
        )
) TO 'wgs_metagenome_accessions.txt' (HEADER FALSE);
