using HarmonyLib;
using UnityEngine;

namespace PatchedConicFixes
{
    // Startup.Instantly + true: load before any scene, persist forever.
    // This is the standard pattern for Harmony-based KSP mods — patches need
    // to be applied as early as possible and must survive scene transitions.
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
