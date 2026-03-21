import os

def get_sequence_from_pdb(filepath):
    """Витягує послідовність амінокислот (3-літерні коди) з PDB файлу за C-альфа атомами."""
    sequence = []
    try:
        with open(filepath, 'r') as f:
            for line in f:
                # Шукаємо рядки атомів, де назва атома - CA (C-alpha)
                if line.startswith("ATOM") and line[12:16].strip() == "CA":
                    res_name = line[17:20].strip()
                    sequence.append(res_name)
    except Exception as e:
        print(f"Помилка читання {filepath}: {e}")
    return sequence

def compare_pdb_folders(alpha_dir, beta_dir):
    # Перевіряємо, чи існують вказані папки
    if not os.path.exists(alpha_dir) or not os.path.exists(beta_dir):
        print("Помилка: Переконайтеся, що папки 'alpha' та 'beta' існують за вказаним шляхом.")
        return

    # Отримуємо список всіх файлів з суфіксом _beta.pdb у папці beta
    beta_files = [f for f in os.listdir(beta_dir) if f.endswith("_beta.pdb")]
    
    if not beta_files:
        print(f"У папці '{beta_dir}' не знайдено файлів з розширенням _beta.pdb")
        return

    print(f"Знайдено {len(beta_files)} файлів у папці '{beta_dir}'. Починаю аналіз...\n")
    print("-" * 60)

    for beta_file in beta_files:
        # Формуємо очікувану назву відповідного alpha файлу
        # Наприклад, з "n_res0_beta.pdb" робимо "n_res0_alpha.pdb"
        base_name = beta_file.replace("_beta.pdb", "")
        alpha_file = base_name + "_alpha.pdb"
        
        beta_path = os.path.join(beta_dir, beta_file)
        alpha_path = os.path.join(alpha_dir, alpha_file)
        
        # 1. Перевіряємо наявність відповідного файлу в папці alpha
        if not os.path.exists(alpha_path):
            print(f"❌ ВІДСУТНІЙ: Для '{beta_file}' немає файлу '{alpha_file}' у папці '{alpha_dir}'")
            continue
            
        # 2. Якщо пара є, зчитуємо та порівнюємо їхні послідовності
        seq_beta = get_sequence_from_pdb(beta_path)
        seq_alpha = get_sequence_from_pdb(alpha_path)
        
        if not seq_beta or not seq_alpha:
            print(f"❓ ПОМИЛКА: Не вдалося прочитати координати CA для пари '{base_name}'")
            continue
            
        if seq_beta == seq_alpha:
            print(f"✅ ЗБІГ: '{base_name}' (Пара знайдена, послідовності ідентичні)")
        else:
            print(f"⚠️ УВАГА: Для '{base_name}' пара знайдена, але ПОСЛІДОВНОСТІ ВІДРІЗНЯЮТЬСЯ!")

    print("-" * 60)
    print("Перевірку завершено.")

# === НАЛАШТУВАННЯ ШЛЯХІВ ===
# Вкажіть тут шляхи до ваших папок. 
# Якщо скрипт лежить у тій самій директорії, де знаходяться папки 'alpha' та 'beta', 
# залиште просто їхні назви:
DIR_ALPHA = "alpha" 
DIR_BETA = "beta"

# Запуск головної функції
if __name__ == "__main__":
    compare_pdb_folders(DIR_ALPHA, DIR_BETA)