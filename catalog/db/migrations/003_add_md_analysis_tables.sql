-- システム情報テーブル (MD計算の系情報を保存)
CREATE TABLE IF NOT EXISTS md_systems (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    vina_result_id INTEGER NOT NULL,
    
    -- システム構成
    total_atoms INTEGER,
    protein_atoms INTEGER,
    ligand_atoms INTEGER,
    water_molecules INTEGER,
    ions_count INTEGER,
    
    -- Box情報
    box_x REAL,
    box_y REAL,
    box_z REAL,
    box_volume REAL,
    
    -- 力場情報
    force_field TEXT,
    water_model TEXT,
    
    -- 計算条件
    temperature REAL DEFAULT 300.0,
    pressure REAL DEFAULT 1.0,
    timestep_ps REAL DEFAULT 0.002,
    
    -- 性能指標
    performance_nsday REAL,
    gpu_model TEXT,
    
    FOREIGN KEY (vina_result_id) REFERENCES vina_results(rowid)
);

-- MD解析結果テーブル
CREATE TABLE IF NOT EXISTS md_analysis (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    vina_result_id INTEGER NOT NULL,
    
    -- RMSD
    rmsd_backbone_avg REAL,
    rmsd_backbone_std REAL,
    rmsd_ligand_avg REAL,
    rmsd_ligand_std REAL,
    
    -- 水素結合
    hbond_avg REAL,
    hbond_std REAL,
    hbond_max INTEGER,
    
    -- Radius of gyration
    rg_avg REAL,
    rg_std REAL,
    
    -- RMSF
    rmsf_avg REAL,
    rmsf_max REAL,
    
    -- エネルギー統計
    potential_energy_avg REAL,
    potential_energy_std REAL,
    
    -- 解析ファイルパス
    analysis_json TEXT,
    
    FOREIGN KEY (vina_result_id) REFERENCES vina_results(rowid)
);

CREATE INDEX IF NOT EXISTS idx_md_systems_vina ON md_systems(vina_result_id);
CREATE INDEX IF NOT EXISTS idx_md_analysis_vina ON md_analysis(vina_result_id);
