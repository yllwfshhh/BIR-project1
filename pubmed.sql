-- Create a table
CREATE TABLE IF NOT EXISTS users (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    name TEXT,
    age INTEGER,
    email TEXT
);

-- Insert data
INSERT INTO users (name, age, email) VALUES ('Alice', 25, 'alice@example.com');

-- Query data
SELECT * FROM users;

-- Update data
UPDATE users SET age = 26 WHERE name = 'Alice';

-- Delete data
DELETE FROM users WHERE name = 'Alice';
