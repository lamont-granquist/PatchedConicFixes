using HarmonyLib;
using UnityEngine;

namespace PatchedConicFixes
{
    // Minimal Harmony proof-of-concept: postfix on PatchedConics._CalculatePatch.
    // The delegate field PatchedConics.CalculatePatch points here by default;
    // _CalculatePatch is the actual static method Harmony can patch.
    [HarmonyPatch(typeof(PatchedConics), nameof(PatchedConics._CalculatePatch))]
    class PatchedConics__CalculatePatch
    {
        static void Postfix(Orbit p, bool __result)
        {
            Debug.Log($"[PatchedConicFixes] _CalculatePatch: result={__result} referenceBody={p?.referenceBody?.name}");
        }
    }
}