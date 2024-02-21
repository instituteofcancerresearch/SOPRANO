from datetime import datetime

from SOPRANO.utils.mpi_utils import as_single_process


def time_output():
    now = datetime.now()
    return now.strftime("%d/%m/%Y %H:%M:%S")


def task_output(msg):
    print(f"[{time_output()}] {msg}", flush=True)


def line_output(n=60):
    print("-" * n, flush=True)


def title_output():
    print(
        """
███████  ██████  ██████  ██████   █████  ███    ██  ██████  
██      ██    ██ ██   ██ ██   ██ ██   ██ ████   ██ ██    ██ 
███████ ██    ██ ██████  ██████  ███████ ██ ██  ██ ██    ██ 
     ██ ██    ██ ██      ██   ██ ██   ██ ██  ██ ██ ██    ██ 
███████  ██████  ██      ██   ██ ██   ██ ██   ████  ██████  
""",
        flush=True,
    )
    line_output()
    print("Selection On PRotein ANnotated regiOns", flush=True)
    line_output()


@as_single_process()
def startup_output(**kwargs):
    title_output()

    if kwargs:
        # print params to console
        for k, v in kwargs.items():
            print("-> {0:.<30}".format(k) + f"{v}", flush=True)

        line_output()
