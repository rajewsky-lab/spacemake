class SpacemakeError(Exception):
    def __init__(self, msg=None):
        self.msg = msg

    def __str__(self):
        msg = 'ERROR: ' + str(self.__class__.__name__) + '\n'

        if hasattr(self, 'msg') and self.msg is not None:
            msg += self.msg

        return msg

class FileWrongExtensionError(SpacemakeError):
    def __init__(self, filename, expected_extension):
        self.filename = filename
        self.expected_extension = expected_extension

    def __str__(self):
        msg = super().__str__()
        msg += f'File {self.filename} has wrong extension.\n'
        msg += f'The extension should be {self.expected_extension}.\n'

        return msg

class ConfigVariableError(SpacemakeError):
    def __init__(self, variable_name, variable_value):
        self.variable_name = variable_name
        self.variable_value = variable_value

class UnrecognisedConfigVariable(SpacemakeError):
    def __init__(self, variable_name, variable_options):
        self.variable_name = variable_name
        self.variable_options = variable_options

    def __str__(self):
        msg = super().__str__()
        msg += f'unrecognised variable {self.variable_name}\n'
        msg += f'it has to be one of {self.variable_options}.'

        return msg

class EmptyConfigVariableError(SpacemakeError):
    def __init__(self, variable_name):
        self.variable_name = variable_name

    def __str__(self):
        msg = super().__str__()
        msg += f'cannot remove, or set {self.variable_name} to emtpy list, or None\n'
        msg += 'this ERROR could happen in two cases: \n'
        msg += f'1) you tried to remove a {self.variable_name}, '
        msg += f'and as a result the sample would not have'
        msg += f' any {self.variable_name} available.\n'
        msg += f'2) you tried to remove the `default` value of'
        msg += f' {self.variable_name} from the configuration.\n'

        return msg

class ConfigVariableNotFoundError(ConfigVariableError):
    def __str__(self):
        msg = super().__str__()
        msg += f'{self.variable_name}: {self.variable_value} not found.\n'
        msg += f'you can add a new {self.variable_name} using the '
        msg += f'`spacemake config add_{self.variable_name}` command.\n'

        return msg

class ConfigVariableIncompleteError(ConfigVariableError):
    def __init__(self, missing_key, **kwargs):
        super().__init__(**kwargs)
        self.missing_key = missing_key

    def __str__(self):
        msg = super().__str__()
        msg += f'{self.variable_name}: {self.variable_value} '
        msg += f'is missing required key {self.required_key}.\n'
        msg += f'You can update this key of {self.variable_value} using the '
        msg += f'`spacemake config update_{self.variable_name}` command.\n'

        return msg

class InvalidBarcodeStructureError(SpacemakeError):
    def __init__(self, tag_name, to_match):
        self.tag_name = tag_name
        self.to_match = to_match

    def __str__(self):
        msg = super().__str__()
        msg += f'{self.tag_name} does not match {self.to_match}.\n'
        msg += f'Example matching would be: r1[0:12] for the first 12n of Read1 '
        msg += f'for {self.tag_name}\n'


class DuplicateConfigVariableError(ConfigVariableError):
    def __str__(self):
        msg = super().__str__()
        msg += f'{self.variable_name}: {self.variable_value} already exists.\n'
        msg += f'To update it use `spacemake config update_{self.variable_name}`,\n'
        msg += f'To delete it use `spacemake config delete_{self.variable_name}.\n'

        return msg

class NoProjectSampleProvidedError(SpacemakeError):
    def __init__(self):
        pass

    def __str__(self):
        msg = super().__str__()
        msg += f'no projects or samples were provided.\n'

        return msg

class ProjectSampleNotFoundError(SpacemakeError):
    def __init__(self, var_name, var_value):
        self.var_name = var_name
        self.var_value = var_value

    def __str__(self):
        msg = super().__str__()
        msg += f'sample with {self.var_name}={self.var_value} not found.\n'
        msg += 'you can add a new sample with `spacemake projects add_sample` command.\n'

        return msg

class SampleAlreadyExistsError(SpacemakeError):
    def __init__(self, ix):
        self.ix = ix

    def __str__(self):
        msg = super().__str__()
        msg += f'sample with (project_id, sample_id)={self.ix} already exists.\n'
        msg += 'in order to update this sample use `spacemake projects update_sample`,\n'
        msg += 'to delete it use `spacemake projects delete_sample`.\n'

        return msg

class InconsistentVariablesDuringMerge(ConfigVariableError):
    def __init__(self, ix, **kwargs):
        super().__init__(**kwargs)
        self.ix = ix

    def __str__(self):
        msg = super().__str__()
        msg += f'\nthe samples that you trying to merge have different '
        msg += f'{self.variable_name} values.\n\ninconsistent values:'
        msg += f' {self.variable_value}\n'
        msg += f'samples: {self.ix}.\n\n'
        msg += 'You can only merge samples which have the same '
        msg += f'{self.variable_name}, or if there is an overlap.\n'

        return msg
