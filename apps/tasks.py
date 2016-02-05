import logging
from multiprocessing import Process, JoinableQueue, Queue, cpu_count

logging.basicConfig()
_logger = logging.getLogger(__name__)


class Poison:  # Sentinel to signal queue end
    pass


class Worker (Process):

    def __init__(self, tasks, results):
        super(Worker, self).__init__()
        self.tasks = tasks
        self.results = results

    def run(self):
        proc_name = self.name
        while True:
            next_task = self.tasks.get()
            if isinstance(next_task, Poison):
                _logger.debug("{}: Exiting".format(proc_name))
                self.tasks.task_done()
                break
            _logger.debug("{}: {}".format(proc_name, next_task))
            res = next_task()
            self.tasks.task_done()
            self.results.put(res)


class Workers (list):
    """
    Runs multiple funcion calls on the SAME computer, using multiprocessing.
    """
    def __init__(self, tasks, results, cpus):
        for x in xrange(cpus):
            self.append(Worker(tasks, results))

    def start(self):
        for pi in self:
            pi.start()


class Task (object):
    """
    A functor object that passes in the function and arguments.
    """
    def __init__(self, target, args):
        self.target = target
        self.args = listify(args)

    def __call__(self):
        return self.target(*self.args)


class Tasks (object):
    """
    Runs multiple function calls, but collects return values.

    Producer-consumer model. See also:
    https://pymotw.com/2/multiprocessing/communication.html
    """
    def __init__(self, target, args, cpus=cpu_count()):
        tasks = JoinableQueue()
        results = Queue()
        njobs = len(args)
        cpus = min(cpus, njobs)

        self.workers = Workers(tasks, results, cpus)
        self.workers.start()
        # Dumping data
        for a in args:
            tasks.put(Task(target, a))
        # Sentry to signal termination
        for i in xrange(cpus):
            tasks.put(Poison())

        # Wait for all tasks to finish
        tasks.join()

        self.results = []
        while njobs:
            res = results.get()
            if isinstance(res, Poison):
                continue
            self.results.append(res)
            njobs -= 1


def listify(a):
    return a if (isinstance(a, list) or isinstance(a, tuple)) else [a]
