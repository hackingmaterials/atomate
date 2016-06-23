from tqdm import tqdm

__author__ = 'Anubhav Jain <ajain@lbl.gov>'


class FileMaterialsBuilder:
    def __init__(self, materials_write, data_file, delimiter=",", header_lines=0):
        """
        Updates the database using a data file. Format of file must be:
        <material_id>, <property>, <value>
        for which <property> is the materials key to update.

        Comment lines should start with '#'.

        Args:
            materials_write: mongodb collection for materials (write access needed)
            data_file: (str) path to data file
            **kwargs: **kwargs for csv reader
        """
        self._materials = materials_write
        self._data_file = data_file
        self._delimiter = delimiter
        self.header_lines = header_lines


    def run(self):
        print("Starting FileMaterials Builder.")
        with open(self._data_file, 'rb') as f:
            line_no = 0
            for line in tqdm(f):
                line = line.strip()
                if not line.startswith("#"):
                    line_no += 1
                    if line_no > self.header_lines:
                        print("Processing line: {}".format(line))
                        line = line.split(self._delimiter)
                        m_id = int(line[0])
                        key = line[1]
                        val = line[2]
                        try:
                            val = float(val)
                        except:
                            pass

                        x = self._materials.find_one_and_update(
                            {"material_id": m_id}, {"$set": {key: val}},
                            {"material_id": 1})
                        if not x:
                            raise ValueError("Could not find "
                                             "material with material_id: {}".
                                             format(m_id))

        print("FileMaterials Builder finished processing")