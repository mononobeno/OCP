-- OCP results database schema

PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS targets (
  id            INTEGER PRIMARY KEY AUTOINCREMENT,
  code          TEXT NOT NULL UNIQUE,    -- "aps_apoh" など
  name          TEXT NOT NULL,           -- 表示名
  pdb_id        TEXT,                    -- 1C1Z, 6C2W ...
  receptor_path TEXT                     -- catalog/targets/.../receptor.pdb
);

CREATE TABLE IF NOT EXISTS libraries (
  id            INTEGER PRIMARY KEY AUTOINCREMENT,
  code          TEXT NOT NULL UNIQUE,    -- "zinc_2d_smi_v1" など
  description   TEXT
);

CREATE TABLE IF NOT EXISTS ligands (
  id            INTEGER PRIMARY KEY AUTOINCREMENT,
  zinc_id       TEXT NOT NULL,
  library_id    INTEGER NOT NULL,
  smiles_hash   TEXT,
  UNIQUE (zinc_id, library_id),
  FOREIGN KEY (library_id) REFERENCES libraries(id)
);

CREATE TABLE IF NOT EXISTS runs (
  id            INTEGER PRIMARY KEY AUTOINCREMENT,
  run_uuid      TEXT NOT NULL UNIQUE,
  started_at    TEXT NOT NULL,
  target_id     INTEGER NOT NULL,
  library_id    INTEGER NOT NULL,
  notes         TEXT,
  FOREIGN KEY (target_id)  REFERENCES targets(id),
  FOREIGN KEY (library_id) REFERENCES libraries(id)
);

CREATE TABLE IF NOT EXISTS vina_results (
  run_id        INTEGER NOT NULL,
  ligand_id     INTEGER NOT NULL,
  mode_rank     INTEGER NOT NULL,
  affinity_kcal REAL,
  rmsd_lb       REAL,
  rmsd_ub       REAL,
  out_relpath   TEXT,
  PRIMARY KEY (run_id, ligand_id, mode_rank),
  FOREIGN KEY (run_id)    REFERENCES runs(id),
  FOREIGN KEY (ligand_id) REFERENCES ligands(id)
);

CREATE TABLE IF NOT EXISTS md_results (
  run_id         INTEGER NOT NULL,
  ligand_id      INTEGER NOT NULL,
  stage          TEXT NOT NULL,   -- "em", "eq", "prod" など
  rmsd_bb        REAL,
  energy_min     REAL,
  stable_flag    INTEGER,
  report_relpath TEXT,
  PRIMARY KEY (run_id, ligand_id, stage),
  FOREIGN KEY (run_id)  REFERENCES runs(id),
  FOREIGN KEY (ligand_id) REFERENCES ligands(id)
);

CREATE TABLE IF NOT EXISTS pages_publish (
  id            INTEGER PRIMARY KEY AUTOINCREMENT,
  run_id        INTEGER NOT NULL,
  published_at  TEXT NOT NULL,
  commit_hash   TEXT,
  url           TEXT,
  FOREIGN KEY (run_id) REFERENCES runs(id)
);
