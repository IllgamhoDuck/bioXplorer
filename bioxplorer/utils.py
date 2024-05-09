

def update_code(file_path, old_line, new_line):
    """Updates a specific line in a file with a new line."""
    with open(file_path, 'r') as file:
        lines = file.readlines()
    lines = [line.replace(old_line, new_line) if line.strip() == old_line else line for line in lines]
    with open(file_path, 'w') as file:
        file.writelines(lines)
