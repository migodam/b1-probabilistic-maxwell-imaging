# Green base-LU completion phase schema: independent read-only review

Snapshot UTC: 2026-10-23T17:40:48.629419+00:00. Only the root ledger and audit source were read; no Book, physics or main was called. No original ledger, source, budget or earlier FAIL snapshot was modified. This is execution-accounting evidence, not scientific GO.

## Result

```json
{
  "snapshot_rows": 72807,
  "snapshot_last_utc": "2026-10-23T17:40:47Z",
  "Green_rows": 1300,
  "reserve_actions": 296,
  "completed_actions": 295,
  "completed_rhs": 1532,
  "missing_phase_completions": 295,
  "ambiguous_or_unreserved_completions": 0,
  "unconsumed_base_action_reserves": 1,
  "duplicate_reserve_ids": 0,
  "max_base_action_residual": 2.066317862862781e-12
}
```

The ordered replay uses the exact stage/PID/scope/label role, an earlier still-unconsumed base_LU_reused reserve, compatible frequency when supplied, and exact RHS quantity. A missing phase is tolerated ONLY for base_LU_action_residual. Each completion must have exactly one eligible reservation and one phase-bearing context; the reservation is then consumed once. No nearest-neighbor guess, future reserve, pooled-context merge or synthetic payment is permitted. This test is stronger than unique-context matching because it also requires a unique individual reservation.

A repeated Krylov label with fresh reserve IDs is not by itself a duplicate payment: each of these actual base-LU calls has its own earlier reserve and its own single completion. Reserve-ID multiplicity was checked against the whole root ledger, not just this scope. Literal re-execution cannot be excluded from ledger alone, but no duplicate charge/receipt or unreserved completed action appears in this ordered accounting evidence.

The scope is serial; if later live rows create an ambiguity or a reserve without completion, those rows must remain ambiguous/in-progress rather than inherit this snapshot conclusion. No claim is made about unsaved/unlogged computations.

The snapshot has **one still-unconsumed live base-action reservation**. It is not a completed action and is not evidence of an unreserved solve; its payment remains retained. Its final outcome requires a later snapshot. The 295 completed actions are all uniquely paired; no ambiguity, missing prior payment, repeated reserve ID or multiply consumed reservation was found among them.

## Per-role totals

| Role | Completed actions | Reserved/completed RHS | Missing completion phase | First reserve line | Last completion line |
|---|---:|---:|---:|---:|---:|
|ROM_primal_R0_4source|1|4|1|71542|71544|
|ROM_volume_goal_dual_12receiver|1|12|1|71546|71548|
|shared_R0V_current_and_enriched_field_ROM|1|16|1|71550|71552|
|EXACT_WOODBURY_R0U_chunk64|20|1228|20|71559|71618|
|KRYLOV_base_LU_apply_1RHS|272|272|272|71650|72805|

## Exceptions / remaining reservations

```json
{
  "anomalies": [],
  "remaining": [
    {
      "reserve_line": 72806,
      "remaining_rhs": 1,
      "label": "KRYLOV_base_LU_apply_1RHS"
    }
  ],
  "duplicate_green_reserve_ids": []
}
```

## Complete ordered pairing evidence

Root-ledger line numbers are one-based physical text lines. Every listed reserve ID is globally unique in the snapshot. Completion phase is missing in the rows explicitly labelled below; inherited phase is context reconstruction for the audit only, never a ledger edit.

| Reserve line | Completion line | Role | RHS | Unique earlier reserve ID | Reserve phase | Completion phase | Intervening same-role operation |
|---:|---:|---|---:|---|---|---|---:|
|71542|71544|ROM_primal_R0_4source|4|85d80907d0704a0a8bca6856632f6f5f|P1B_green_comparison|MISSING|1|
|71546|71548|ROM_volume_goal_dual_12receiver|12|2e32d25451a547d99d1bd73f33d4e319|P1B_green_comparison|MISSING|1|
|71550|71552|shared_R0V_current_and_enriched_field_ROM|16|44565d50188f42d29b42451efc0ed85a|P1B_green_comparison|MISSING|1|
|71559|71561|EXACT_WOODBURY_R0U_chunk64|64|d5faae7f32d44d1e856b56e22b5fa8ee|P1B_green_comparison|MISSING|1|
|71562|71564|EXACT_WOODBURY_R0U_chunk64|64|85e6c13ef57f4fa5b888361950b2608d|P1B_green_comparison|MISSING|1|
|71565|71567|EXACT_WOODBURY_R0U_chunk64|64|a3da1d1dbba54a9da7290ae11a084069|P1B_green_comparison|MISSING|1|
|71568|71570|EXACT_WOODBURY_R0U_chunk64|64|e210f3e7ef2c4ba79d82b6a23e25e6b6|P1B_green_comparison|MISSING|1|
|71571|71573|EXACT_WOODBURY_R0U_chunk64|64|8a7adb3d501640f3839a9553b2fac0d0|P1B_green_comparison|MISSING|1|
|71574|71576|EXACT_WOODBURY_R0U_chunk64|64|5e4a8a1e64a44af0a20c820544912a7f|P1B_green_comparison|MISSING|1|
|71577|71579|EXACT_WOODBURY_R0U_chunk64|64|8cb8f7363d61473d87259e6181db0c32|P1B_green_comparison|MISSING|1|
|71580|71582|EXACT_WOODBURY_R0U_chunk64|64|df3d27d30a5541e98692595042dc5575|P1B_green_comparison|MISSING|1|
|71583|71585|EXACT_WOODBURY_R0U_chunk64|64|87af31bdf9c74cb7a7e85fb23cb6beb4|P1B_green_comparison|MISSING|1|
|71586|71588|EXACT_WOODBURY_R0U_chunk64|64|d489c51f13584cb79f33d5927ec6605f|P1B_green_comparison|MISSING|1|
|71589|71591|EXACT_WOODBURY_R0U_chunk64|64|f00ab2efe1074c3b851a2ddf4a3f59e9|P1B_green_comparison|MISSING|1|
|71592|71594|EXACT_WOODBURY_R0U_chunk64|64|91e0e5307dcd42918d67385fe17670af|P1B_green_comparison|MISSING|1|
|71595|71597|EXACT_WOODBURY_R0U_chunk64|64|635dc034521848b2a40a0256faa47dd4|P1B_green_comparison|MISSING|1|
|71598|71600|EXACT_WOODBURY_R0U_chunk64|64|17cddfb5e4a64e8c87c1033ddff6caf5|P1B_green_comparison|MISSING|1|
|71601|71603|EXACT_WOODBURY_R0U_chunk64|64|79cc8d669e404114843a5d17b7435f52|P1B_green_comparison|MISSING|1|
|71604|71606|EXACT_WOODBURY_R0U_chunk64|64|2bf546c375594a73b7e34dfd320ed056|P1B_green_comparison|MISSING|1|
|71607|71609|EXACT_WOODBURY_R0U_chunk64|64|52958dcd2b3b4b9aac22f644efbac283|P1B_green_comparison|MISSING|1|
|71610|71612|EXACT_WOODBURY_R0U_chunk64|64|8e1127598f254c689246676594c44210|P1B_green_comparison|MISSING|1|
|71613|71615|EXACT_WOODBURY_R0U_chunk64|64|481861a4301641348372d6f4b5b1d2bf|P1B_green_comparison|MISSING|1|
|71616|71618|EXACT_WOODBURY_R0U_chunk64|12|5234ce6b2b0149e5a4d4acf206bbd09e|P1B_green_comparison|MISSING|1|
|71650|71652|KRYLOV_base_LU_apply_1RHS|1|3ee8491561314b60bb08c150c25d36fb|P1B_green_comparison|MISSING|1|
|71653|71655|KRYLOV_base_LU_apply_1RHS|1|9a01865cdbec4b90884a5135c6626bae|P1B_green_comparison|MISSING|1|
|71656|71658|KRYLOV_base_LU_apply_1RHS|1|825fe474fe7c423abc2b3f8bd1ce3113|P1B_green_comparison|MISSING|1|
|71659|71661|KRYLOV_base_LU_apply_1RHS|1|fc9d5df99aa043149d5580711b9ba825|P1B_green_comparison|MISSING|1|
|71662|71664|KRYLOV_base_LU_apply_1RHS|1|cb24c55b62ca433f8063f2617f58345a|P1B_green_comparison|MISSING|1|
|71665|71667|KRYLOV_base_LU_apply_1RHS|1|c6d168ab635647aca3f965b89626bf47|P1B_green_comparison|MISSING|1|
|71668|71670|KRYLOV_base_LU_apply_1RHS|1|d330ffe98f2c44d795c55dbf1f83df5b|P1B_green_comparison|MISSING|1|
|71671|71673|KRYLOV_base_LU_apply_1RHS|1|1fcf4e02493046509d733370792d0134|P1B_green_comparison|MISSING|1|
|71674|71676|KRYLOV_base_LU_apply_1RHS|1|88bf50332063414e92ac9bce4c13e29a|P1B_green_comparison|MISSING|1|
|71677|71679|KRYLOV_base_LU_apply_1RHS|1|a5a24e1e60d64327a67faf65328bebcc|P1B_green_comparison|MISSING|1|
|71680|71682|KRYLOV_base_LU_apply_1RHS|1|8b94e475ed0d431e93bb2555f41e8530|P1B_green_comparison|MISSING|1|
|71683|71685|KRYLOV_base_LU_apply_1RHS|1|c0a0949f161b4a77878340a1da3461da|P1B_green_comparison|MISSING|1|
|71686|71688|KRYLOV_base_LU_apply_1RHS|1|de0ba91004574f80ad3151bd110d9d60|P1B_green_comparison|MISSING|1|
|71689|71691|KRYLOV_base_LU_apply_1RHS|1|c33847ad5b1c4674b804544b00f81b8d|P1B_green_comparison|MISSING|1|
|71692|71694|KRYLOV_base_LU_apply_1RHS|1|46f6bd2a066d4f3b894195da6901df21|P1B_green_comparison|MISSING|1|
|71695|71697|KRYLOV_base_LU_apply_1RHS|1|5b80f286f87e48beac697d4c356b3cc7|P1B_green_comparison|MISSING|1|
|71698|71700|KRYLOV_base_LU_apply_1RHS|1|1d035adb68c14375b79524386e789ae9|P1B_green_comparison|MISSING|1|
|71701|71703|KRYLOV_base_LU_apply_1RHS|1|37d862b4ad904f56a2a7ba6389edeaed|P1B_green_comparison|MISSING|1|
|71704|71706|KRYLOV_base_LU_apply_1RHS|1|be175059a55647f9a69cff911b8a7f40|P1B_green_comparison|MISSING|1|
|71707|71709|KRYLOV_base_LU_apply_1RHS|1|f29e8154f14f4ccc99c91e9c5dea789a|P1B_green_comparison|MISSING|1|
|71710|71712|KRYLOV_base_LU_apply_1RHS|1|d7461e27107347e9b2fc3ec18d5209d0|P1B_green_comparison|MISSING|1|
|71713|71715|KRYLOV_base_LU_apply_1RHS|1|5f80c96ecded4d04a88d48f5f0d8152b|P1B_green_comparison|MISSING|1|
|71716|71718|KRYLOV_base_LU_apply_1RHS|1|2e382633e18b48b2a5f7076299809c81|P1B_green_comparison|MISSING|1|
|71719|71721|KRYLOV_base_LU_apply_1RHS|1|92e1d9c2c82e45de84571e943aa96ab2|P1B_green_comparison|MISSING|1|
|71722|71724|KRYLOV_base_LU_apply_1RHS|1|700a99d6258c4393ab259255aba3d0d4|P1B_green_comparison|MISSING|1|
|71725|71727|KRYLOV_base_LU_apply_1RHS|1|d838d42c7eb643fc95a9596087a2fc22|P1B_green_comparison|MISSING|1|
|71728|71730|KRYLOV_base_LU_apply_1RHS|1|c2da5f984d8341f5a7d144d01c1d7ce8|P1B_green_comparison|MISSING|1|
|71731|71733|KRYLOV_base_LU_apply_1RHS|1|94d5974008c34474ae210b16653b4026|P1B_green_comparison|MISSING|1|
|71734|71736|KRYLOV_base_LU_apply_1RHS|1|aa1edc97a91745cb8e29763d06620ca5|P1B_green_comparison|MISSING|1|
|71737|71739|KRYLOV_base_LU_apply_1RHS|1|62e482a617db4e7c9ff05c1affb5632f|P1B_green_comparison|MISSING|1|
|71740|71742|KRYLOV_base_LU_apply_1RHS|1|9d005be946e24188906dfb89d8194e1d|P1B_green_comparison|MISSING|1|
|71743|71745|KRYLOV_base_LU_apply_1RHS|1|0a4312575e9e4950a228b67a5b79b57a|P1B_green_comparison|MISSING|1|
|71780|71782|KRYLOV_base_LU_apply_1RHS|1|c36f62ec8a3d4e598246dcc806a0972c|P1B_green_comparison|MISSING|1|
|71783|71785|KRYLOV_base_LU_apply_1RHS|1|3af622235e8d4882bd5cbee71b888fe0|P1B_green_comparison|MISSING|1|
|71786|71788|KRYLOV_base_LU_apply_1RHS|1|439b9e9a99fb40e29340028116b7169e|P1B_green_comparison|MISSING|1|
|71789|71791|KRYLOV_base_LU_apply_1RHS|1|2b0ee39c3b6f4550b73307e370619004|P1B_green_comparison|MISSING|1|
|71792|71794|KRYLOV_base_LU_apply_1RHS|1|804346c8d77948e5b3a5588a859260f2|P1B_green_comparison|MISSING|1|
|71795|71797|KRYLOV_base_LU_apply_1RHS|1|55bb351657a2411183a7d9d01359da1e|P1B_green_comparison|MISSING|1|
|71798|71800|KRYLOV_base_LU_apply_1RHS|1|3f36d0f4a0224bdb8934de12e3f521e8|P1B_green_comparison|MISSING|1|
|71801|71803|KRYLOV_base_LU_apply_1RHS|1|7cb6ba69a02d434687ba7fe0fbf27698|P1B_green_comparison|MISSING|1|
|71804|71806|KRYLOV_base_LU_apply_1RHS|1|4ae0fdfd403c48c68fe8b4d7c2a139bd|P1B_green_comparison|MISSING|1|
|71807|71809|KRYLOV_base_LU_apply_1RHS|1|9ad4cc18aa8f4489b71ad5ec9c5bc2f6|P1B_green_comparison|MISSING|1|
|71810|71812|KRYLOV_base_LU_apply_1RHS|1|0b8d5a7460ac4f7790e3589b17075117|P1B_green_comparison|MISSING|1|
|71813|71815|KRYLOV_base_LU_apply_1RHS|1|6251711d7ffc4ad49dbfdd8f30ccfa0e|P1B_green_comparison|MISSING|1|
|71816|71818|KRYLOV_base_LU_apply_1RHS|1|dc30b3c607ed4caf8cea9c2102bf43f5|P1B_green_comparison|MISSING|1|
|71819|71821|KRYLOV_base_LU_apply_1RHS|1|0d1d7f2d2f644e6aa1427f3dab3b3496|P1B_green_comparison|MISSING|1|
|71822|71824|KRYLOV_base_LU_apply_1RHS|1|a91a916fe4014d9ea3e723953113fa62|P1B_green_comparison|MISSING|1|
|71825|71827|KRYLOV_base_LU_apply_1RHS|1|c90907e3d3b745c095f121c7bb468007|P1B_green_comparison|MISSING|1|
|71828|71830|KRYLOV_base_LU_apply_1RHS|1|862733da65af4f929b3db55cb9bcf3a3|P1B_green_comparison|MISSING|1|
|71831|71833|KRYLOV_base_LU_apply_1RHS|1|cbda62970e014f44b82d681033792be5|P1B_green_comparison|MISSING|1|
|71834|71836|KRYLOV_base_LU_apply_1RHS|1|64b62c3af8ed46a7aeda6c9f19e20535|P1B_green_comparison|MISSING|1|
|71837|71839|KRYLOV_base_LU_apply_1RHS|1|e3bc6c1d07574ceda89c5414c087e893|P1B_green_comparison|MISSING|1|
|71840|71842|KRYLOV_base_LU_apply_1RHS|1|9cae4627cdf346519046a35ed7734db4|P1B_green_comparison|MISSING|1|
|71843|71845|KRYLOV_base_LU_apply_1RHS|1|69c09c5de4fb41cdac41053080e6ba8c|P1B_green_comparison|MISSING|1|
|71846|71848|KRYLOV_base_LU_apply_1RHS|1|c86ba87efdf944fba18c64cd5c961fdb|P1B_green_comparison|MISSING|1|
|71849|71851|KRYLOV_base_LU_apply_1RHS|1|5e1b2b5cff3142a1b0d6ba89a295826a|P1B_green_comparison|MISSING|1|
|71852|71854|KRYLOV_base_LU_apply_1RHS|1|7dbcfe83dd4141da9fad2f0cdfd4bdbf|P1B_green_comparison|MISSING|1|
|71855|71857|KRYLOV_base_LU_apply_1RHS|1|d732237322ef48f6b193d8ee119ffd9f|P1B_green_comparison|MISSING|1|
|71858|71860|KRYLOV_base_LU_apply_1RHS|1|8a7a98a4d34b4de4b97c903fe76d7c59|P1B_green_comparison|MISSING|1|
|71861|71863|KRYLOV_base_LU_apply_1RHS|1|c71992d6d7c842919a2b1e45dcb737a6|P1B_green_comparison|MISSING|1|
|71864|71866|KRYLOV_base_LU_apply_1RHS|1|4331e948891948059d6189442ab52756|P1B_green_comparison|MISSING|1|
|71867|71869|KRYLOV_base_LU_apply_1RHS|1|5a6f64beb95e49d4931a9ad9e034575f|P1B_green_comparison|MISSING|1|
|71870|71872|KRYLOV_base_LU_apply_1RHS|1|334315fb381c44418dd1fe8b974dcdb3|P1B_green_comparison|MISSING|1|
|71873|71875|KRYLOV_base_LU_apply_1RHS|1|150d38aebab24be6a2176a7d65284f03|P1B_green_comparison|MISSING|1|
|71910|71912|KRYLOV_base_LU_apply_1RHS|1|9f5d4475182b4852993d6bba7648ebfe|P1B_green_comparison|MISSING|1|
|71913|71915|KRYLOV_base_LU_apply_1RHS|1|bad2b6d7de6543dbbe032221bac11511|P1B_green_comparison|MISSING|1|
|71916|71918|KRYLOV_base_LU_apply_1RHS|1|a4681dc3773b48aa843564c1488cfcaf|P1B_green_comparison|MISSING|1|
|71919|71921|KRYLOV_base_LU_apply_1RHS|1|a6032f28f49b4452922453d36f121de3|P1B_green_comparison|MISSING|1|
|71922|71924|KRYLOV_base_LU_apply_1RHS|1|c6e46f9ac5834e5ebce088af2f2652bd|P1B_green_comparison|MISSING|1|
|71925|71927|KRYLOV_base_LU_apply_1RHS|1|b436ad1436bf4ff68c015d9a89d63367|P1B_green_comparison|MISSING|1|
|71928|71930|KRYLOV_base_LU_apply_1RHS|1|a384b85e11b44c5a8a94a0d574a95a7a|P1B_green_comparison|MISSING|1|
|71931|71933|KRYLOV_base_LU_apply_1RHS|1|fb0ef7f563dd4c0caad3b0ba1bba2e23|P1B_green_comparison|MISSING|1|
|71934|71936|KRYLOV_base_LU_apply_1RHS|1|e9eaade9d7eb4e4b8856cc63be598941|P1B_green_comparison|MISSING|1|
|71937|71939|KRYLOV_base_LU_apply_1RHS|1|cad8d9f141d0473c9c1cd070270c226a|P1B_green_comparison|MISSING|1|
|71940|71942|KRYLOV_base_LU_apply_1RHS|1|f330310ca0a642aea0a03d2021deabc7|P1B_green_comparison|MISSING|1|
|71943|71945|KRYLOV_base_LU_apply_1RHS|1|66834a265c7b4257a87130fed27d0bc5|P1B_green_comparison|MISSING|1|
|71946|71948|KRYLOV_base_LU_apply_1RHS|1|c209f8b1bfcb4f2ebf742883b519f8a0|P1B_green_comparison|MISSING|1|
|71949|71951|KRYLOV_base_LU_apply_1RHS|1|ae76a89f07a7472492ee87b519ad2ce1|P1B_green_comparison|MISSING|1|
|71952|71954|KRYLOV_base_LU_apply_1RHS|1|f5e95ad88f0e4bbe90fcb727b7192d00|P1B_green_comparison|MISSING|1|
|71955|71957|KRYLOV_base_LU_apply_1RHS|1|6c99e34b9b814e62a1b96d3076eb41f9|P1B_green_comparison|MISSING|1|
|71958|71960|KRYLOV_base_LU_apply_1RHS|1|cbbd910f5ab24bc982d49e24f28fc5df|P1B_green_comparison|MISSING|1|
|71961|71963|KRYLOV_base_LU_apply_1RHS|1|dffea7cf657a4a8b8876d30311f1fe9b|P1B_green_comparison|MISSING|1|
|71964|71966|KRYLOV_base_LU_apply_1RHS|1|ce4e15ef32d1481cb6fc62734ead869c|P1B_green_comparison|MISSING|1|
|71967|71969|KRYLOV_base_LU_apply_1RHS|1|2be32bfe4f454f7ab4ff6a38c3da4d98|P1B_green_comparison|MISSING|1|
|71970|71972|KRYLOV_base_LU_apply_1RHS|1|837a1724bf40463fa086bd1241f81e14|P1B_green_comparison|MISSING|1|
|71973|71975|KRYLOV_base_LU_apply_1RHS|1|12509129aa004125aa9be028b3002e1d|P1B_green_comparison|MISSING|1|
|71976|71978|KRYLOV_base_LU_apply_1RHS|1|136048da70064c92b795f0dfed0267fc|P1B_green_comparison|MISSING|1|
|71979|71981|KRYLOV_base_LU_apply_1RHS|1|a4f919a5452c49dc83556891f714ce93|P1B_green_comparison|MISSING|1|
|72016|72018|KRYLOV_base_LU_apply_1RHS|1|25cc7efc6fda405e9a6ba91371b14895|P1B_green_comparison|MISSING|1|
|72019|72021|KRYLOV_base_LU_apply_1RHS|1|988a6d0e9fb345b3b64188c5992706b4|P1B_green_comparison|MISSING|1|
|72022|72024|KRYLOV_base_LU_apply_1RHS|1|10d3d20868be4e259cb44911db2167ad|P1B_green_comparison|MISSING|1|
|72025|72027|KRYLOV_base_LU_apply_1RHS|1|6b460e4e1341497eb813c96fbcb91b76|P1B_green_comparison|MISSING|1|
|72028|72030|KRYLOV_base_LU_apply_1RHS|1|827da9f36a024294a977f521680880f7|P1B_green_comparison|MISSING|1|
|72031|72033|KRYLOV_base_LU_apply_1RHS|1|575ede2a74f44da9a588770dd5f9d569|P1B_green_comparison|MISSING|1|
|72034|72036|KRYLOV_base_LU_apply_1RHS|1|be739dc5d3dd411b821bc8e22afb05c6|P1B_green_comparison|MISSING|1|
|72037|72039|KRYLOV_base_LU_apply_1RHS|1|744995e628bc480692e1fd582ab92682|P1B_green_comparison|MISSING|1|
|72040|72042|KRYLOV_base_LU_apply_1RHS|1|42a8e30c30c845399dd441903e17d309|P1B_green_comparison|MISSING|1|
|72043|72045|KRYLOV_base_LU_apply_1RHS|1|3ccd49ac5a9c4fffb9a57b67c9d4e253|P1B_green_comparison|MISSING|1|
|72046|72048|KRYLOV_base_LU_apply_1RHS|1|5fbee1193e634f6d9bbe456573a7dae5|P1B_green_comparison|MISSING|1|
|72049|72051|KRYLOV_base_LU_apply_1RHS|1|9360ae1e0f6c4f8eb3566092044a3422|P1B_green_comparison|MISSING|1|
|72052|72054|KRYLOV_base_LU_apply_1RHS|1|647811eeff1748acb141f80a19a2cf38|P1B_green_comparison|MISSING|1|
|72055|72057|KRYLOV_base_LU_apply_1RHS|1|deb3ad7349394a299ab5d375f29bc6ca|P1B_green_comparison|MISSING|1|
|72058|72060|KRYLOV_base_LU_apply_1RHS|1|a9ef363e015841089304f989e56a5901|P1B_green_comparison|MISSING|1|
|72061|72063|KRYLOV_base_LU_apply_1RHS|1|b7e0ba9195794a2097b6b956e11ab714|P1B_green_comparison|MISSING|1|
|72064|72066|KRYLOV_base_LU_apply_1RHS|1|f81bd2ec9a164109a00fd4b6e8345d8f|P1B_green_comparison|MISSING|1|
|72067|72069|KRYLOV_base_LU_apply_1RHS|1|0e66a01935704ce4a837c9eeb4c90544|P1B_green_comparison|MISSING|1|
|72070|72072|KRYLOV_base_LU_apply_1RHS|1|e1b24b95c77745aba59fa3b32dc5248d|P1B_green_comparison|MISSING|1|
|72073|72075|KRYLOV_base_LU_apply_1RHS|1|af846e0dffe947f3ade0e425fc74315f|P1B_green_comparison|MISSING|1|
|72076|72078|KRYLOV_base_LU_apply_1RHS|1|e8b0a4ef02074d81bcd9e322373ac80d|P1B_green_comparison|MISSING|1|
|72079|72081|KRYLOV_base_LU_apply_1RHS|1|3ab1eb13f8f349f9b8a7b4a7a76603cf|P1B_green_comparison|MISSING|1|
|72082|72084|KRYLOV_base_LU_apply_1RHS|1|c495139bd7f540769956e66edd0fd9ab|P1B_green_comparison|MISSING|1|
|72085|72087|KRYLOV_base_LU_apply_1RHS|1|2649d1b969114142a5525fcb3b989b6c|P1B_green_comparison|MISSING|1|
|72122|72124|KRYLOV_base_LU_apply_1RHS|1|598e4940771f48ff87a0f68a383f9335|P1B_green_comparison|MISSING|1|
|72125|72127|KRYLOV_base_LU_apply_1RHS|1|d6676cc32a83416ea074da6860aa8df2|P1B_green_comparison|MISSING|1|
|72128|72130|KRYLOV_base_LU_apply_1RHS|1|0e94ccbceda146fc8b4ba8ccefd281b6|P1B_green_comparison|MISSING|1|
|72131|72133|KRYLOV_base_LU_apply_1RHS|1|d58988429a0142168457fee6ff3027b6|P1B_green_comparison|MISSING|1|
|72134|72136|KRYLOV_base_LU_apply_1RHS|1|42f5f84a3e9144bf8408ab08be3a2c50|P1B_green_comparison|MISSING|1|
|72137|72139|KRYLOV_base_LU_apply_1RHS|1|2d9873165e4c4a7eb1c64b87fb35f95a|P1B_green_comparison|MISSING|1|
|72140|72142|KRYLOV_base_LU_apply_1RHS|1|ba22dff94eba4bb6b6bfc671d5ac28e7|P1B_green_comparison|MISSING|1|
|72143|72145|KRYLOV_base_LU_apply_1RHS|1|c3c7bb7887a2413f9a5ce6b2ee962498|P1B_green_comparison|MISSING|1|
|72146|72148|KRYLOV_base_LU_apply_1RHS|1|54a68dd048ee4676a6bbcefbdad8c539|P1B_green_comparison|MISSING|1|
|72149|72151|KRYLOV_base_LU_apply_1RHS|1|522bd870405c49a59800a2ed25d5490c|P1B_green_comparison|MISSING|1|
|72152|72154|KRYLOV_base_LU_apply_1RHS|1|cd4b1402df9c4358a83ed81b49e0e072|P1B_green_comparison|MISSING|1|
|72155|72157|KRYLOV_base_LU_apply_1RHS|1|c2c0489194504a13b13ffe66c9993ec1|P1B_green_comparison|MISSING|1|
|72158|72160|KRYLOV_base_LU_apply_1RHS|1|2b8b3353e1a544be96a2b64bef60cef1|P1B_green_comparison|MISSING|1|
|72161|72163|KRYLOV_base_LU_apply_1RHS|1|f44293de21454b63b25e8eab0892e3d4|P1B_green_comparison|MISSING|1|
|72164|72166|KRYLOV_base_LU_apply_1RHS|1|743115d304964d44955a33bf3fb77a8e|P1B_green_comparison|MISSING|1|
|72167|72169|KRYLOV_base_LU_apply_1RHS|1|ab4c20f6bf32427d828eaf077ae8d0b5|P1B_green_comparison|MISSING|1|
|72170|72172|KRYLOV_base_LU_apply_1RHS|1|2d036cf55d9a45bbb59497b8a4fb6f73|P1B_green_comparison|MISSING|1|
|72173|72175|KRYLOV_base_LU_apply_1RHS|1|c0deb2c3028d41b1b408536d84b178bd|P1B_green_comparison|MISSING|1|
|72176|72178|KRYLOV_base_LU_apply_1RHS|1|7284e968f57449ccb8965bd23805f01a|P1B_green_comparison|MISSING|1|
|72179|72181|KRYLOV_base_LU_apply_1RHS|1|d6f26c8355e64032b66b5c0da1932c51|P1B_green_comparison|MISSING|1|
|72182|72184|KRYLOV_base_LU_apply_1RHS|1|40785748b2b64ff6a2661537c7901bd6|P1B_green_comparison|MISSING|1|
|72185|72187|KRYLOV_base_LU_apply_1RHS|1|b810c40b15394d818c2871f2f068c948|P1B_green_comparison|MISSING|1|
|72188|72190|KRYLOV_base_LU_apply_1RHS|1|94179ff913c7475cb51187eb66fe9859|P1B_green_comparison|MISSING|1|
|72191|72193|KRYLOV_base_LU_apply_1RHS|1|a8cf307b1b254c2e92f45a4bbef15505|P1B_green_comparison|MISSING|1|
|72228|72230|KRYLOV_base_LU_apply_1RHS|1|2c5ce764b3c449f8b6df3615fa86aa52|P1B_green_comparison|MISSING|1|
|72231|72233|KRYLOV_base_LU_apply_1RHS|1|6903ff2d724d4a2b8057d43b0c5261c9|P1B_green_comparison|MISSING|1|
|72234|72236|KRYLOV_base_LU_apply_1RHS|1|d127b39327f4450eb55df3b52c059acc|P1B_green_comparison|MISSING|1|
|72237|72239|KRYLOV_base_LU_apply_1RHS|1|aed0a5aede444e70a16e33505f497ca1|P1B_green_comparison|MISSING|1|
|72240|72242|KRYLOV_base_LU_apply_1RHS|1|0b4e9d0187eb46279b623b681501dac2|P1B_green_comparison|MISSING|1|
|72243|72245|KRYLOV_base_LU_apply_1RHS|1|37f1bb93c45143b593faba80e5573972|P1B_green_comparison|MISSING|1|
|72246|72248|KRYLOV_base_LU_apply_1RHS|1|47b8b395b2094f459e101be47683a241|P1B_green_comparison|MISSING|1|
|72249|72251|KRYLOV_base_LU_apply_1RHS|1|e2f9eada77484d37b20a0a8eca70f185|P1B_green_comparison|MISSING|1|
|72252|72254|KRYLOV_base_LU_apply_1RHS|1|dc98f0ec7abb426f8086add649dd2623|P1B_green_comparison|MISSING|1|
|72255|72257|KRYLOV_base_LU_apply_1RHS|1|43fd64f775ce4d9f97ae4cc1b81236bb|P1B_green_comparison|MISSING|1|
|72258|72260|KRYLOV_base_LU_apply_1RHS|1|5148137c61c34e92b9448a4475a28c2e|P1B_green_comparison|MISSING|1|
|72261|72263|KRYLOV_base_LU_apply_1RHS|1|e6e2fc55339b4220b9fb552d3b30d6ea|P1B_green_comparison|MISSING|1|
|72264|72266|KRYLOV_base_LU_apply_1RHS|1|29667e79f0cd4feab19f0a13c1c16824|P1B_green_comparison|MISSING|1|
|72267|72269|KRYLOV_base_LU_apply_1RHS|1|7e54bb4e768047fdbd74867fe8e1c6d7|P1B_green_comparison|MISSING|1|
|72270|72272|KRYLOV_base_LU_apply_1RHS|1|0f96b9ef40a64150bbeb483f4633b8d7|P1B_green_comparison|MISSING|1|
|72273|72275|KRYLOV_base_LU_apply_1RHS|1|3185e0c24d1d4f10a359c9ac55300ba8|P1B_green_comparison|MISSING|1|
|72276|72278|KRYLOV_base_LU_apply_1RHS|1|82d9172df04544a991a354756a7e2269|P1B_green_comparison|MISSING|1|
|72279|72281|KRYLOV_base_LU_apply_1RHS|1|4dc0783682484cf59e672ffa1c406ef3|P1B_green_comparison|MISSING|1|
|72282|72284|KRYLOV_base_LU_apply_1RHS|1|ab7b6a837b6947d18e3c02fbffb2cf58|P1B_green_comparison|MISSING|1|
|72285|72287|KRYLOV_base_LU_apply_1RHS|1|fd30b0a995a8498aa81b73423cdf0e21|P1B_green_comparison|MISSING|1|
|72288|72290|KRYLOV_base_LU_apply_1RHS|1|9571a6b2ace242fdb7422d17297d4f18|P1B_green_comparison|MISSING|1|
|72291|72293|KRYLOV_base_LU_apply_1RHS|1|928136130cbc4c7a8e42b77609baa705|P1B_green_comparison|MISSING|1|
|72294|72296|KRYLOV_base_LU_apply_1RHS|1|b0f7dceebe0b4b0e938c59cc78d47327|P1B_green_comparison|MISSING|1|
|72297|72299|KRYLOV_base_LU_apply_1RHS|1|4c3d4b35acca461585dc90ca21f9c3e8|P1B_green_comparison|MISSING|1|
|72334|72336|KRYLOV_base_LU_apply_1RHS|1|64f242a9ba9744298805cab820dad73e|P1B_green_comparison|MISSING|1|
|72337|72339|KRYLOV_base_LU_apply_1RHS|1|c1768ec9777a4d8e9f478592dd97708a|P1B_green_comparison|MISSING|1|
|72340|72342|KRYLOV_base_LU_apply_1RHS|1|37e9f360a8a641d2a58b01a1b550221b|P1B_green_comparison|MISSING|1|
|72343|72345|KRYLOV_base_LU_apply_1RHS|1|3bd104a20d9a45f59a1678c5178e6f98|P1B_green_comparison|MISSING|1|
|72346|72348|KRYLOV_base_LU_apply_1RHS|1|73d7d3013e4446d6a88598c030f4e650|P1B_green_comparison|MISSING|1|
|72349|72351|KRYLOV_base_LU_apply_1RHS|1|a33c1acd7a4c43d4b844add5e5206fe7|P1B_green_comparison|MISSING|1|
|72352|72354|KRYLOV_base_LU_apply_1RHS|1|23ffba29a88d43aca63a0a8090893812|P1B_green_comparison|MISSING|1|
|72355|72357|KRYLOV_base_LU_apply_1RHS|1|19a0ae439be045adaffb4225a9c3b35d|P1B_green_comparison|MISSING|1|
|72358|72360|KRYLOV_base_LU_apply_1RHS|1|3eb0167c5a4c4f70bd24a0d8e4d9b6e5|P1B_green_comparison|MISSING|1|
|72361|72363|KRYLOV_base_LU_apply_1RHS|1|a22caf27c8304471bfef1170ccc2b344|P1B_green_comparison|MISSING|1|
|72364|72366|KRYLOV_base_LU_apply_1RHS|1|1c3a333184ba42cebd86424805dbeb7f|P1B_green_comparison|MISSING|1|
|72367|72369|KRYLOV_base_LU_apply_1RHS|1|263f32daf60b44e7aa7903fb81e75ec2|P1B_green_comparison|MISSING|1|
|72370|72372|KRYLOV_base_LU_apply_1RHS|1|140172e809eb423e85c45b570b9d42bc|P1B_green_comparison|MISSING|1|
|72373|72375|KRYLOV_base_LU_apply_1RHS|1|0c50a1bfda174450afe14be3e236084a|P1B_green_comparison|MISSING|1|
|72376|72378|KRYLOV_base_LU_apply_1RHS|1|69499124710849a2bfc9de887c53a854|P1B_green_comparison|MISSING|1|
|72379|72381|KRYLOV_base_LU_apply_1RHS|1|ca13e3d4de024087804a4d0c406d174e|P1B_green_comparison|MISSING|1|
|72382|72384|KRYLOV_base_LU_apply_1RHS|1|22ffc7b85b4b46ff911944e3397e7dfa|P1B_green_comparison|MISSING|1|
|72385|72387|KRYLOV_base_LU_apply_1RHS|1|8bf2aba2d6584373adb482380602a6cb|P1B_green_comparison|MISSING|1|
|72388|72390|KRYLOV_base_LU_apply_1RHS|1|c6320595354e4ad1b2442fa8ce237fb1|P1B_green_comparison|MISSING|1|
|72391|72393|KRYLOV_base_LU_apply_1RHS|1|4093b2cf6a8446059cd46a3861c477bc|P1B_green_comparison|MISSING|1|
|72394|72396|KRYLOV_base_LU_apply_1RHS|1|cd1090c7e1304327b6b78c0adbb2915a|P1B_green_comparison|MISSING|1|
|72397|72399|KRYLOV_base_LU_apply_1RHS|1|d4e1d53c96894b079160ef33f4ecdfe2|P1B_green_comparison|MISSING|1|
|72400|72402|KRYLOV_base_LU_apply_1RHS|1|4aa21941b28e47d5b3123c56daa94f97|P1B_green_comparison|MISSING|1|
|72403|72405|KRYLOV_base_LU_apply_1RHS|1|630c9a972f3041c7bd2e9dc3c8b311a0|P1B_green_comparison|MISSING|1|
|72406|72408|KRYLOV_base_LU_apply_1RHS|1|b859260932584f618ddb92612a114447|P1B_green_comparison|MISSING|1|
|72409|72411|KRYLOV_base_LU_apply_1RHS|1|736ba106108a48538e18aeb642873c44|P1B_green_comparison|MISSING|1|
|72412|72414|KRYLOV_base_LU_apply_1RHS|1|051c27fd010b4cf2b93ca63239520bfd|P1B_green_comparison|MISSING|1|
|72415|72417|KRYLOV_base_LU_apply_1RHS|1|3aabbdcb1f4c4de2894a846ed90d7c6c|P1B_green_comparison|MISSING|1|
|72452|72454|KRYLOV_base_LU_apply_1RHS|1|26ae1b98284b482aa3f4e445d41a8269|P1B_green_comparison|MISSING|1|
|72455|72457|KRYLOV_base_LU_apply_1RHS|1|914378df669b4b99abe61dccd317b7ce|P1B_green_comparison|MISSING|1|
|72458|72460|KRYLOV_base_LU_apply_1RHS|1|c1a90548bb16475eaca8516442bd800e|P1B_green_comparison|MISSING|1|
|72461|72463|KRYLOV_base_LU_apply_1RHS|1|e6e865a56efc40b7ba496dff18105622|P1B_green_comparison|MISSING|1|
|72464|72466|KRYLOV_base_LU_apply_1RHS|1|a9c3166ed50a48d6b698f5e6039680a5|P1B_green_comparison|MISSING|1|
|72467|72469|KRYLOV_base_LU_apply_1RHS|1|5962790e96d24148b4a4486d0c93b5de|P1B_green_comparison|MISSING|1|
|72470|72472|KRYLOV_base_LU_apply_1RHS|1|f97727e6bd244d19a332c39a4c7f473a|P1B_green_comparison|MISSING|1|
|72473|72475|KRYLOV_base_LU_apply_1RHS|1|0bb4feee261d4d7fb042d66456a6b745|P1B_green_comparison|MISSING|1|
|72476|72478|KRYLOV_base_LU_apply_1RHS|1|899ef8e9b87a4758ae065d0fe7674721|P1B_green_comparison|MISSING|1|
|72479|72481|KRYLOV_base_LU_apply_1RHS|1|c9a1a649910e4e03a47107c90e0764d1|P1B_green_comparison|MISSING|1|
|72482|72484|KRYLOV_base_LU_apply_1RHS|1|f499570a68884ed9af2b47f8e47673a9|P1B_green_comparison|MISSING|1|
|72485|72487|KRYLOV_base_LU_apply_1RHS|1|28148cdd8cda4a80bf9ca6880d11b2b7|P1B_green_comparison|MISSING|1|
|72488|72490|KRYLOV_base_LU_apply_1RHS|1|1c0676e8f9ac4ddab5da408b66a37530|P1B_green_comparison|MISSING|1|
|72491|72493|KRYLOV_base_LU_apply_1RHS|1|8eab60707c824f8f81a1400cac4cf445|P1B_green_comparison|MISSING|1|
|72494|72496|KRYLOV_base_LU_apply_1RHS|1|fa3afe88031244c49470bea955ac6b26|P1B_green_comparison|MISSING|1|
|72497|72499|KRYLOV_base_LU_apply_1RHS|1|0ee78180210343469702413bfa3179f7|P1B_green_comparison|MISSING|1|
|72500|72502|KRYLOV_base_LU_apply_1RHS|1|ce7e53dd06d448ff8a6f35f3882bcbbd|P1B_green_comparison|MISSING|1|
|72503|72505|KRYLOV_base_LU_apply_1RHS|1|3605c4b2825f4367be63368a9c7d6d7d|P1B_green_comparison|MISSING|1|
|72506|72508|KRYLOV_base_LU_apply_1RHS|1|598c0bf8d38f406a99b2513a5bf66424|P1B_green_comparison|MISSING|1|
|72509|72511|KRYLOV_base_LU_apply_1RHS|1|d5dcecf20f3a4aa6aab54b348b5b9247|P1B_green_comparison|MISSING|1|
|72546|72548|KRYLOV_base_LU_apply_1RHS|1|c13f68265786436eac017461ce1e8741|P1B_green_comparison|MISSING|1|
|72549|72551|KRYLOV_base_LU_apply_1RHS|1|54634e404c414850a777363a75976594|P1B_green_comparison|MISSING|1|
|72552|72554|KRYLOV_base_LU_apply_1RHS|1|87f4f185ade144ec9d67c6c43c1f1c63|P1B_green_comparison|MISSING|1|
|72555|72557|KRYLOV_base_LU_apply_1RHS|1|046bbb7b652f431f9643c3e781d52b3e|P1B_green_comparison|MISSING|1|
|72558|72560|KRYLOV_base_LU_apply_1RHS|1|df3ed27ba54e41c6bebca3aa0d195213|P1B_green_comparison|MISSING|1|
|72561|72563|KRYLOV_base_LU_apply_1RHS|1|132a3b60f22c4a0f9b0623e0ab7ba305|P1B_green_comparison|MISSING|1|
|72564|72566|KRYLOV_base_LU_apply_1RHS|1|4a37e29cf2dc4bd1ba4c47e2c7ee68e1|P1B_green_comparison|MISSING|1|
|72567|72569|KRYLOV_base_LU_apply_1RHS|1|2b95978a6c9e4408a5d1e0706daad241|P1B_green_comparison|MISSING|1|
|72570|72572|KRYLOV_base_LU_apply_1RHS|1|22917810e81f44ccbef5c244b492052e|P1B_green_comparison|MISSING|1|
|72573|72575|KRYLOV_base_LU_apply_1RHS|1|ea47bde1a97a44019b041dc57d4fa734|P1B_green_comparison|MISSING|1|
|72576|72578|KRYLOV_base_LU_apply_1RHS|1|a3280eaf5d64419d92eac946b4706fdb|P1B_green_comparison|MISSING|1|
|72579|72581|KRYLOV_base_LU_apply_1RHS|1|ff79f3193ae249148baf10572a573355|P1B_green_comparison|MISSING|1|
|72582|72584|KRYLOV_base_LU_apply_1RHS|1|8809d6a4cc9140d0894e2109f0fba461|P1B_green_comparison|MISSING|1|
|72585|72587|KRYLOV_base_LU_apply_1RHS|1|2b76d27576ec41cebb1aabef7108183d|P1B_green_comparison|MISSING|1|
|72588|72590|KRYLOV_base_LU_apply_1RHS|1|783223b78d3044c49f05206242bc4c6c|P1B_green_comparison|MISSING|1|
|72591|72593|KRYLOV_base_LU_apply_1RHS|1|239a80fbda90467fa99ae97a37ed71d5|P1B_green_comparison|MISSING|1|
|72594|72596|KRYLOV_base_LU_apply_1RHS|1|923ca4e0d4e94ba283345fe39df51bae|P1B_green_comparison|MISSING|1|
|72597|72599|KRYLOV_base_LU_apply_1RHS|1|b06d854ed80148b2ae3368fbb61d64ac|P1B_green_comparison|MISSING|1|
|72600|72602|KRYLOV_base_LU_apply_1RHS|1|b7eb75df911d461489c698f4275209ef|P1B_green_comparison|MISSING|1|
|72603|72605|KRYLOV_base_LU_apply_1RHS|1|27abfdd9c27a407da36be4afd14cd218|P1B_green_comparison|MISSING|1|
|72640|72642|KRYLOV_base_LU_apply_1RHS|1|5e972cd53fc8407c8d8cedaa291b5a49|P1B_green_comparison|MISSING|1|
|72643|72645|KRYLOV_base_LU_apply_1RHS|1|1c0d1a0952014c3dbf575ff786e917e0|P1B_green_comparison|MISSING|1|
|72646|72648|KRYLOV_base_LU_apply_1RHS|1|2bede86ec9d246f387d9d2f29282e9a3|P1B_green_comparison|MISSING|1|
|72649|72651|KRYLOV_base_LU_apply_1RHS|1|9ec870ee7e2547b4bfa1be6232c477c3|P1B_green_comparison|MISSING|1|
|72652|72654|KRYLOV_base_LU_apply_1RHS|1|89748bb7c1ba42048a1fcdd03161b828|P1B_green_comparison|MISSING|1|
|72655|72657|KRYLOV_base_LU_apply_1RHS|1|07accb16e7ee4c79ac29244cf724064d|P1B_green_comparison|MISSING|1|
|72658|72660|KRYLOV_base_LU_apply_1RHS|1|2d4c244c143a47e196679fc76504c09f|P1B_green_comparison|MISSING|1|
|72661|72663|KRYLOV_base_LU_apply_1RHS|1|a63a1b4551ee41969fa6088d5cf0b83e|P1B_green_comparison|MISSING|1|
|72664|72666|KRYLOV_base_LU_apply_1RHS|1|fbc4b8dfaa624625bcf1259eac8f6625|P1B_green_comparison|MISSING|1|
|72667|72669|KRYLOV_base_LU_apply_1RHS|1|a974989a8c3c4389b3d526a630c4d208|P1B_green_comparison|MISSING|1|
|72670|72672|KRYLOV_base_LU_apply_1RHS|1|dc5bbe4f375d464f9b0cc851f65d6135|P1B_green_comparison|MISSING|1|
|72673|72675|KRYLOV_base_LU_apply_1RHS|1|86e261d25a994d8b88722ece6b99b0e4|P1B_green_comparison|MISSING|1|
|72676|72678|KRYLOV_base_LU_apply_1RHS|1|7e35b5b6ed0d4fb6bf9f1f41e997dbf7|P1B_green_comparison|MISSING|1|
|72679|72681|KRYLOV_base_LU_apply_1RHS|1|549716df5e4141b9be876e6a7a239241|P1B_green_comparison|MISSING|1|
|72682|72684|KRYLOV_base_LU_apply_1RHS|1|6b797b3380f84e159f6d3a36e8b151e9|P1B_green_comparison|MISSING|1|
|72685|72687|KRYLOV_base_LU_apply_1RHS|1|f138c22c93ec40dda152525f48e516ae|P1B_green_comparison|MISSING|1|
|72688|72690|KRYLOV_base_LU_apply_1RHS|1|79211f9612bb4321b136f0e9deedf11f|P1B_green_comparison|MISSING|1|
|72691|72693|KRYLOV_base_LU_apply_1RHS|1|4bd0c773221c41d5ab95f4cadef3521b|P1B_green_comparison|MISSING|1|
|72694|72696|KRYLOV_base_LU_apply_1RHS|1|bcffaf8410b94a9f9cc5243555f392f8|P1B_green_comparison|MISSING|1|
|72697|72699|KRYLOV_base_LU_apply_1RHS|1|67851783d7ba42ca8daf3cd4da130cb5|P1B_green_comparison|MISSING|1|
|72700|72702|KRYLOV_base_LU_apply_1RHS|1|6fe53d4c88ca47869d50f14ca4c5afe0|P1B_green_comparison|MISSING|1|
|72703|72705|KRYLOV_base_LU_apply_1RHS|1|f65d44c7acb141de87733260f6fd41bf|P1B_green_comparison|MISSING|1|
|72706|72708|KRYLOV_base_LU_apply_1RHS|1|5c471bb52b8a41e59174c88f92b3ae4a|P1B_green_comparison|MISSING|1|
|72709|72711|KRYLOV_base_LU_apply_1RHS|1|58388189fd8b42a9a3536b4cf1b47795|P1B_green_comparison|MISSING|1|
|72712|72714|KRYLOV_base_LU_apply_1RHS|1|626fe5df8fc94381b41e85d9bf9329a9|P1B_green_comparison|MISSING|1|
|72715|72717|KRYLOV_base_LU_apply_1RHS|1|f5467d1054de4028bac19067063ebe53|P1B_green_comparison|MISSING|1|
|72718|72720|KRYLOV_base_LU_apply_1RHS|1|5aa1e31b41d742ae8d63a2c5d24b50ec|P1B_green_comparison|MISSING|1|
|72721|72723|KRYLOV_base_LU_apply_1RHS|1|ab15d9e431df4278a7d32b334e165e83|P1B_green_comparison|MISSING|1|
|72758|72760|KRYLOV_base_LU_apply_1RHS|1|ff44479a879f4edb93c0e477d9e325e3|P1B_green_comparison|MISSING|1|
|72761|72763|KRYLOV_base_LU_apply_1RHS|1|2c794d57ed634c08915af500d18e2c0a|P1B_green_comparison|MISSING|1|
|72764|72766|KRYLOV_base_LU_apply_1RHS|1|65605b047c1e4ae194e6f42c84eae70f|P1B_green_comparison|MISSING|1|
|72767|72769|KRYLOV_base_LU_apply_1RHS|1|6e64c18ad82b45f180c6b6e9c9f030fe|P1B_green_comparison|MISSING|1|
|72770|72772|KRYLOV_base_LU_apply_1RHS|1|5ac5024bca574034a0414f8cad9ab8ce|P1B_green_comparison|MISSING|1|
|72773|72775|KRYLOV_base_LU_apply_1RHS|1|a451f25baf4749dda7d80513ee349074|P1B_green_comparison|MISSING|1|
|72776|72778|KRYLOV_base_LU_apply_1RHS|1|2ba117d471cf45f388f3b45801170cdb|P1B_green_comparison|MISSING|1|
|72779|72781|KRYLOV_base_LU_apply_1RHS|1|ed0a4b5fd4e8406f91fc09b3b89765a9|P1B_green_comparison|MISSING|1|
|72782|72784|KRYLOV_base_LU_apply_1RHS|1|8b25005b97f540b6b192a536a793db00|P1B_green_comparison|MISSING|1|
|72785|72787|KRYLOV_base_LU_apply_1RHS|1|395608a51423495b8df990efec575dcb|P1B_green_comparison|MISSING|1|
|72788|72790|KRYLOV_base_LU_apply_1RHS|1|34dda92967294278b2f88f3dadd8693d|P1B_green_comparison|MISSING|1|
|72791|72793|KRYLOV_base_LU_apply_1RHS|1|ccd2df2c376c4eab907f0f8b8f32bb86|P1B_green_comparison|MISSING|1|
|72794|72796|KRYLOV_base_LU_apply_1RHS|1|fb9de23c43ab4bd18454922e87ff983a|P1B_green_comparison|MISSING|1|
|72797|72799|KRYLOV_base_LU_apply_1RHS|1|5f6b1243474b43c7b560ef39dd2f75c1|P1B_green_comparison|MISSING|1|
|72800|72802|KRYLOV_base_LU_apply_1RHS|1|32e765bea7cd4684bf14532fc26acde8|P1B_green_comparison|MISSING|1|
|72803|72805|KRYLOV_base_LU_apply_1RHS|1|5ff1eb164794437d8074619a7be2dc83|P1B_green_comparison|MISSING|1|

## Code cause and narrow correction

`src/wr/green_reuse_control.py::_action` reserves through `p.book.reserve(p.phase, label, rhs=...)`, but emits `base_LU_action_residual` through `p.book.event(...)` without passing phase. The reserve and completion still share stage/PID/scope/label and exact RHS. The examples 71542→71544 (4 RHS) and 71546→71548 (12 RHS) are included in the complete evidence above.

The proposed auditor-only correction is consistent with this snapshot if restricted to this completion event with missing phase and a unique unconsumed earlier context. Ambiguity, quantity overrun, duplicate completion, wrong PID/stage/scope/role or explicitly conflicting phase must still fail. Preserve the original FAIL artifact as a schema-diagnostic snapshot; do not relabel or delete its historical finding. This review does not require or authorize changes to ledger, Book, budget or physical source.
