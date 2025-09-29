from pathlib import Path

if __name__ == "__main__":
    modules = Path(__file__).parent.glob("*py")

    modules_names = {m.stem for m in modules} - {"__init__"}
    print("Available modules are: ", modules_names)
