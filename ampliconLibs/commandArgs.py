class CommandArgs:

    def __init__(self, command_args, boolean_flags=[]):
        """
        do not append a positional argument to a boolean (no value) flag argument,
        CommandArgs can not handle this situation unless you specified boolean_flags;
        boolean_flags is not needed if there is not any positional arguments
        """
        self.boolean_flags = boolean_flags
        self.args, self.kwargs = self.parse_args(command_args)

    def parse_args(self, args):
        args_list = args.strip().split()
        args = []
        kwargs = {}
        key_flag = None
        for n, arg in enumerate(args_list):
            if arg.startswith('-'):
                if key_flag:
                    kwargs[key_flag] = ""
                if arg in self.boolean_flags:
                    kwargs[arg] = ""
                    key_flag = None
                else:
                    key_flag = arg
            else:
                if key_flag:
                    kwargs[key_flag] = arg
                    key_flag = None
                else:
                    args.append(arg)
        if key_flag:
            kwargs[key_flag] = ""
        return args, kwargs

    @classmethod
    def extend(cls, args1, args2, boolean_flags=[], ignore_positional=True):
        ca1 = cls(args1, boolean_flags)
        ca2 = cls(args2, boolean_flags)
        args = ca1.args
        if not ignore_positional:
            for n, a in enumerate(ca2.args):
                args[n] = a
        kwargs = ca1.kwargs
        kwargs.update(ca2.kwargs)
        for k, v in kwargs.items():
            args.append(k)
            args.append(v)
        return " ".join(args)
