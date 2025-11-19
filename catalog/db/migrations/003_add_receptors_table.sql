-- Migration: Add receptors table
-- Date: 2025-11-20
-- Description: Create receptors table to store protein target information

CREATE TABLE IF NOT EXISTS receptors (
  id            INTEGER PRIMARY KEY AUTOINCREMENT,
  name          TEXT NOT NULL UNIQUE,  -- e.g., "aps_prothrombin", "aps_apoh"
  pdb_id        TEXT,                  -- e.g., "6C2W", "1C1Z"
  description   TEXT,
  pdb_path      TEXT NOT NULL,         -- Path to receptor.pdb
  pdbqt_path    TEXT,                  -- Path to receptor.pdbqt (for Vina)
  
  -- Docking box parameters (for AutoDock Vina)
  box_center_x  REAL,
  box_center_y  REAL,
  box_center_z  REAL,
  box_size_x    REAL DEFAULT 20.0,
  box_size_y    REAL DEFAULT 20.0,
  box_size_z    REAL DEFAULT 20.0,
  
  created_at    TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- Create index for faster lookups
CREATE INDEX IF NOT EXISTS idx_receptors_name ON receptors(name);
CREATE INDEX IF NOT EXISTS idx_receptors_pdb_id ON receptors(pdb_id);
