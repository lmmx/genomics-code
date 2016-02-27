SELECT
DISTINCT s1.name AS patch,
         s2.name AS chromosome,
	 a.exc_seq_region_start AS start,
	 a.exc_seq_region_end AS end 
  FROM assembly_exception AS a,
       seq_region AS s1,
       seq_region AS s2

 WHERE a.seq_region_id = s1.seq_region_id
   AND a.exc_seq_region_id = s2.seq_region_id
   AND a.exc_type IN ('PATCH_FIX', 'PATCH_NOVEL');
