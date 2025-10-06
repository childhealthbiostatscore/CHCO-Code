import pandas as pd
import uuid

# Load harmonized dataset
harmonized = pd.read_csv(
    "/Users/shivaniramesh/Library/CloudStorage/OneDrive-UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/harmonized_dataset.csv",
    low_memory=False
)

# Get unique MRNs
mrns = harmonized["mrn"].dropna().unique()


def generate_f_uuid():
    u = uuid.uuid4().hex  # random UUID as hex string
    return "f" + u[1:]    # replacing first character with 'f' for fake :)


mrn_uuid_map = {mrn: generate_f_uuid() for mrn in mrns}

map_df = pd.DataFrame(list(mrn_uuid_map.items()), columns=["mrn", "uuid"])
map_df.to_csv("/Users/shivaniramesh/Library/CloudStorage/OneDrive-UW/Laura Pyle's files - Biostatistics Core Shared Drive/Data Harmonization/Data Clean/mrn_uuid_map.csv", index=False)

print(f"Saved {len(map_df)} MRN → UUID mappings to mrn_uuid_map.csv")
