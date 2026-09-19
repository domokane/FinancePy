from pathlib import Path

# ============================================================================
# FINANCEPY EXAMPLES - Rename
# ============================================================================

root = Path(__file__).parent

for file in root.rglob("*.py"):
    # Don't rename this script itself
    if file == Path(__file__):
        continue

    # Don't add the prefix twice
    if file.name.startswith("example_"):
        continue

    new_file = file.with_name(f"example_{file.name}")

    if new_file.exists():
        print(f"SKIPPED: {new_file} already exists")
        continue

    print(f"{file.name} -> {new_file.name}")
    file.rename(new_file)

print("\nDone.")
input("Press Enter to close...")