import Setup as b
import MafProcessing



for key,value in b.data["pw"].items():
    print(f"=== {key} ===")
    print(f"Gap Rate: {value['GapRate']}")
    print(f"Sub Rate: {value['SubstitutionRate']}\n")