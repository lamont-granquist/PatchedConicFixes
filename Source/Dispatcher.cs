using System;
using System.Collections;
using System.Threading;

namespace PatchedConicFixes
{
    public class Dispatcher
    {
        private readonly Thread _mainThread;
        private readonly object _lockObject = new object();
        private readonly Queue  _actions    = new Queue();

        public Dispatcher(Thread mainThread) => _mainThread = mainThread;

        public void InvokeAsync(Action action)
        {
            if (Thread.CurrentThread == _mainThread)
                action();
            else
                lock (_lockObject)
                    _actions.Enqueue(action);
        }

        /// <summary>
        ///     Call this from the main thread to drain and execute queued actions.
        /// </summary>
        public void ProcessActions()
        {
            lock (_lockObject)
                while (_actions.Count > 0)
                    ((Action)_actions.Dequeue())();
        }
    }
}
