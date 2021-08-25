class FileWrongExtensionError(Exception):
    def __init__(self, filename, expected_extension):
        self.filename = filename
        self.expected_extension = expected_extension

    def __str__(self):
        msg = f'File {self.filename} has wrong extension.\n'
        msg += f'The extension should be {self.expected_extension}\n'

        return msg

class BarcodeFlavorNotFoundError(Exception):
    def __init__(self, barcode_flavor):
        self.barcode_flavor = barcode_flavor

    def __str__(self):
        msg = f'ERROR: barcode_flavor: {self.barcode_flavor} is not specified.\n'
        msg += f'you can add a new run_mode with `spacemake config add_barcode_flavor`\n'
        
        return msg

class SpeciesNotFoundError(Exception):
    def __init__(self, species):
        self.species = species

    def __str__(self):
        msg = f'ERROR: barcode_flavor: {self.species} is not a valid species.\n'
        msg += f'you can add a new species with `spacemake config add_species`\n'

        return msg

class RunModeNotFoundError(Exception):
    def __init__(self, run_mode_name):
        self.run_mode_name = run_mode_name

    def __str__(self):
        msg = f'ERROR: run_mode: {self.run_mode_name} not found.\n'
        msg += 'you can add a new run_mode using the `spacemake config add_run_mode` command.\n'

        return msg

class NoProjectSampleProvidedError(Exception):
    def __init__(self):
        pass

    def __str__(self):
        msg = f'ERROR: no project or samples were provided.\n'

        return msg

