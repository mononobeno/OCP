-- tools/sql/debug_vina_entities.sql
-- Vina 周辺のエンティティ状態をざっと確認するための SQL 集合。
-- 必要に応じて列名やテーブル名は実際の schema に合わせて調整してください。
-- 例:
--   sqlite3 catalog/db/ocp_results.sqlite < tools/sql/debug_vina_entities.sql

.echo ON

.print '=== [1] libraries ==='
SELECT id, code, description
FROM libraries
ORDER BY id
LIMIT 20;

.print ''
.print '=== [2] targets ==='
SELECT id, code, name
FROM targets
ORDER BY id
LIMIT 20;

.print ''
.print '=== [3] RDKit 3D 生成状況 (ligands) ==='
SELECT
  COUNT(*) AS total_ligands,
  SUM(CASE WHEN has_3d = 1 THEN 1 ELSE 0 END) AS ligands_with_3d,
  ROUND(
    100.0 * SUM(CASE WHEN has_3d = 1 THEN 1 ELSE 0 END) / COUNT(*),
    2
  ) AS pct_with_3d
FROM ligands;

.print ''
.print '=== [4] RDKit (rdkit_etkdg) で 3D 生成済み ligands 件数 ==='
SELECT
  COUNT(*) AS ligands_rdkit_etkdg
FROM ligands
WHERE has_3d = 1
  AND conformer_method = 'rdkit_etkdg';

.print ''
.print '=== [5] 最近の runs 一覧 ==='
SELECT
  id,
  run_uuid,
  started_at,
  library_id,
  target_id,
  notes
FROM runs
ORDER BY id DESC
LIMIT 10;

.print ''
.print '=== [6] Vina 結果 (vina_results) の概要 ==='
SELECT
  COUNT(*) AS total_results,
  MIN(affinity_kcal) AS best_score,
  MAX(affinity_kcal) AS worst_score
FROM vina_results
WHERE mode_rank = 1;

.print ''
.print '=== [7] Vina 結果ランキング (上位 20 件) ==='
SELECT
  vr.run_id,
  vr.ligand_id,
  l.zinc_id,
  vr.affinity_kcal,
  vr.out_relpath
FROM vina_results vr
JOIN ligands l ON l.id = vr.ligand_id
WHERE vr.mode_rank = 1
ORDER BY vr.affinity_kcal ASC
LIMIT 20;

.print ''
.print '=== [8] RUN ごとの結果件数 ==='
SELECT
  r.id AS run_id,
  r.notes,
  COUNT(vr.ligand_id) AS result_count,
  MIN(vr.affinity_kcal) AS best_score,
  MAX(vr.affinity_kcal) AS worst_score
FROM runs r
LEFT JOIN vina_results vr ON vr.run_id = r.id AND vr.mode_rank = 1
GROUP BY r.id, r.notes
ORDER BY r.id DESC
LIMIT 20;

.print ''
.print '=== [END] debug_vina_entities.sql finished ==='
