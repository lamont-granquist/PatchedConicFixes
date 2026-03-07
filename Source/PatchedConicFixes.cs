using System.Threading;
using HarmonyLib;
using UnityEngine;

namespace PatchedConicFixes
{
    // Startup.Instantly + true: load before any scene, persist forever.
    [KSPAddon(KSPAddon.Startup.Instantly, true)]
    public class PatchedConicFixes : MonoBehaviour
    {
        private static readonly Dispatcher _dispatcher = new Dispatcher(Thread.CurrentThread);

        private static void SafePrint(object message) => _dispatcher.InvokeAsync(() => Debug.Log("[PatchedConicFixes] " + message));

        static PatchedConicFixes()
        {
            Logger.GlobalRegister(SafePrint);
        }

        /*
        private void Awake()
        {
            Logger.GlobalRegister(SafePrint);
        }
        */

        private void Start()
        {
            var harmony = new Harmony("PatchedConicFixes");
            harmony.PatchAll();
            Debug.Log("[PatchedConicFixes] Harmony patches applied.");
        }

        private void Update()
        {
            _dispatcher.ProcessActions();
        }
    }
}
