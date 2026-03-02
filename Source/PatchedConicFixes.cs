using HarmonyLib;
using UnityEngine;

namespace PatchedConicFixes
{
    // Startup.Instantly + true: load before any scene, persist forever.
    [KSPAddon(KSPAddon.Startup.Instantly, true)]
    public class PatchedConicFixes : MonoBehaviour
    {
        private void Start()
        {
            var harmony = new Harmony("PatchedConicFixes");
            harmony.PatchAll();
            Debug.Log("[PatchedConicFixes] Harmony patches applied.");
        }
    }
}
