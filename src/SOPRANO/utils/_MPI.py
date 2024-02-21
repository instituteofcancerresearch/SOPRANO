class COMM_WORLD:
    @staticmethod
    def Barrier():
        pass

    @staticmethod
    def Get_rank():
        return 0

    @staticmethod
    def Get_size():
        return 1

    @staticmethod
    def send(*args, **kwargs):
        pass

    @staticmethod
    def recv(*args, **kwargs):
        raise RuntimeError("Invalid operation for single mpi process.")


def Finalize():
    pass
